# Paulo Nuin March 2018


def get_barcodes(barcode_file):


    barcodes_dict ={}
    barcodes = open(barcode_file).read().splitlines()

    for barcode in barcodes:
        temp = barcode.split('\t')
        barcodes_dict[temp[0]] = temp[1]

    return barcodes_dict



if __name__ == '__main__':

    pipeline_barcodes = get_barcodes('all_barcodes.txt')
    legacy_barcodes = get_barcodes('excel_barcodes_OA.txt')

    all_ids = set(legacy_barcodes.keys()).union(set(pipeline_barcodes.keys()))

    matches = []
    not_matched = []
    for i in legacy_barcodes:
        try:
            matches.append((i, legacy_barcodes[i], pipeline_barcodes[i]))
        except:
            not_matched.append((i, legacy_barcodes[i], pipeline_barcodes.has_key(i)))

    print('Total barcodes in Excel file:\t\t\t\t' + str(len(legacy_barcodes)))
    print('Total barcodes from pipeline:\t\t\t\t' + str(len(pipeline_barcodes)))

    print('Total matches:\t\t\t\t\t\t\t\t' + str(len(matches)))
    print('Total not matched (ids not in pipeline)\t\t' + str(len(not_matched)) + '\n')

    print('Matches')
    print('Sample ID\tBarcode Excel\tBarcode pipeline')
    for item in matches:
        print(item[0] + '\t' + item[1] + '\t' + item[2])

    print('\n')

    print('Not matched')
    print('Sample ID\tBarcode Excel\tSample ID in pipeline')
    for item in not_matched:
        print(item[0] + '\t' + item[1] + '\t' + str(item[2]))

    # for i in legacy_barcodes:
    #     print(i)

    # index = 0
    # shared_barcodes = set(pipeline_barcodes.items()) & set(legacy_barcodes.items())
    # for i in pipeline_barcodes:
    #     try:  
    #     except:
    #         print(i)

    # print(index)
    # print(len(shared_barcodes))