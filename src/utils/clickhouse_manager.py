import os
import click
from dotenv import load_dotenv
from clickhouse_driver import Client
from rich.console import Console
from rich.table import Table
from datetime import datetime

# Load environment variables
load_dotenv()

# ClickHouse connection details
CH_HOST = os.getenv('DB_HOST', '192.168.1.167')
CH_PORT = int(os.getenv('DB_PORT', 9000))
CH_USER = os.getenv('DB_USER', 'default')
CH_PASSWORD = os.getenv('DB_PASSWORD', '')
CH_DATABASE = os.getenv('DB_NAME', 'cnv_database')

console = Console()


def get_clickhouse_client():
    return Client(host=CH_HOST, port=CH_PORT, user=CH_USER, password=CH_PASSWORD, database=CH_DATABASE)


def print_table(headers, rows):
    table = Table(show_header=True, header_style="bold magenta")
    for header in headers:
        table.add_column(header)
    for row in rows:
        table.add_row(*[str(item) for item in row])
    console.print(table)


@click.group()
def cli():
    """Manage and monitor ClickHouse CNV database."""
    pass


@cli.command()
def setup():
    """Create database and table if they don't exist."""
    client = get_clickhouse_client()

    # Create database if not exists
    client.execute(f"CREATE DATABASE IF NOT EXISTS {CH_DATABASE}")

    # Create table if not exists
    client.execute(f"""
        CREATE TABLE IF NOT EXISTS {CH_DATABASE}.cnv_data (
            sample_id String,
            run_id String,
            gene String,
            location String,
            count UInt32,
            timestamp DateTime
        ) ENGINE = MergeTree()
        ORDER BY (sample_id, run_id, gene)
    """)

    console.print("[green]Database and table created/verified successfully.[/green]")


@cli.command()
@click.option('--sample-id', required=True, help='Sample ID')
@click.option('--run-id', required=True, help='Run ID')
@click.option('--cnv-file', required=True, type=click.Path(exists=True), help='Path to CNV file')
def insert(sample_id, run_id, cnv_file):
    """Insert sample data into the database."""
    client = get_clickhouse_client()

    with open(cnv_file, 'r') as f:
        next(f)  # Skip header
        data = []
        for line in f:
            location, count = line.strip().split('\t')
            gene = location.split(':')[0]
            data.append((sample_id, run_id, gene, location, int(count), datetime.now()))

    client.execute(f"""
        INSERT INTO {CH_DATABASE}.cnv_data 
        (sample_id, run_id, gene, location, count, timestamp) VALUES
    """, data)

    console.print(f"[green]Inserted {len(data)} records for sample {sample_id}[/green]")


@cli.command()
def monitor():
    """Display database statistics."""
    client = get_clickhouse_client()

    # Get total number of rows
    total_rows = client.execute(f"SELECT COUNT(*) FROM {CH_DATABASE}.cnv_data")[0][0]

    # Get number of samples
    num_samples = client.execute(f"SELECT COUNT(DISTINCT sample_id) FROM {CH_DATABASE}.cnv_data")[0][0]

    # Get number of runs
    num_runs = client.execute(f"SELECT COUNT(DISTINCT run_id) FROM {CH_DATABASE}.cnv_data")[0][0]

    # Get size of the table
    table_size = client.execute(
        f"SELECT formatReadableSize(total_bytes) FROM system.tables WHERE database = '{CH_DATABASE}' AND name = 'cnv_data'")[
        0][0]

    # Get recent inserts
    recent_inserts = client.execute(f"""
        SELECT COUNT(*) 
        FROM {CH_DATABASE}.cnv_data 
        WHERE timestamp >= toDateTime(now() - INTERVAL 1 DAY)
    """)[0][0]

    console.print("[bold]Database Statistics:[/bold]")
    print_table(
        ["Metric", "Value"],
        [
            ["Total Rows", total_rows],
            ["Number of Samples", num_samples],
            ["Number of Runs", num_runs],
            ["Table Size", table_size],
            ["Inserts in Last 24h", recent_inserts]
        ]
    )


@cli.command()
def query():
    """Perform sample queries on the database."""
    client = get_clickhouse_client()

    console.print("\n[bold]Sample Queries:[/bold]")

    # Average count per gene
    avg_count_query = f"""
        SELECT gene, AVG(count) as avg_count
        FROM {CH_DATABASE}.cnv_data
        GROUP BY gene
        ORDER BY avg_count DESC
        LIMIT 10
    """
    result = client.execute(avg_count_query)
    console.print("\nTop 10 genes by average count:")
    print_table(["Gene", "Avg Count"], result)

    # Samples with highest total count
    high_count_query = f"""
        SELECT sample_id, SUM(count) as total_count
        FROM {CH_DATABASE}.cnv_data
        GROUP BY sample_id
        ORDER BY total_count DESC
        LIMIT 5
    """
    result = client.execute(high_count_query)
    console.print("\nTop 5 samples by total count:")
    print_table(["Sample ID", "Total Count"], result)


if __name__ == "__main__":
    cli()