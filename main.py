import csv
import multiprocessing
import time
import requests


def normalize_name(name):
    return name.capitalize().replace('_', ' ')


def worker(accession, data, data_failed):
    attributes = []
    bogus_values = ['missing', 'n/a', 'not provided', 'not provided; submitted under MIGS 2.1', '', 'NA', '/', 'N/A',
                    'missed', 'not located', 'Not Applicable', 'Missing', 'unknown', 'Unknown', 'not collected',
                    'Not applicable', 'not applicable']
    detail_request_url = 'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/' + accession + '/dataset_report?filters.assembly_version=all_assemblies&page_size=1000'
    try:
        result_detail = requests.get(url=detail_request_url)
    except:
        data_failed.append(accession)
        return

    # RefSeq
    attributes.append(['RefSeq', result_detail.json()['reports'][0]['accession']])

    # Sample details
    biosample = result_detail.json()['reports'][0]['assembly_info']['biosample']
    if biosample['accession'] not in bogus_values:
        attributes.append(['BioSample ID', biosample['accession']])
    if 'description' in biosample:
        if biosample['description']['title'] not in bogus_values:
            attributes.append(['Description', biosample['description']['title']])
        if 'comment' in biosample['description']:
            if biosample['description']['comment'] not in bogus_values:
                attributes.append(['Comment', biosample['description']['comment']])

    if 'owner' in biosample:
        if biosample['owner']['name'] not in bogus_values:
            attributes.append(['Owner name', biosample['owner']['name']])

    for item in biosample['attributes']:
        if 'value' in item:
            if item['value'] not in bogus_values:
                attributes.append([normalize_name(item['name']), item['value']])

    if 'models' in biosample:
        attributes.append(['Models', ', '.join(biosample['models'])])

    if 'sample_ids' in biosample:
        for item in biosample['sample_ids']:
            if 'label' in item:
                if item['value'] not in bogus_values:
                    attributes.append([normalize_name(item['label']), item['value']])
            elif 'db' in item:
                if item['value'] not in bogus_values:
                    attributes.append([normalize_name(item['db']), item['value']])

    if biosample['package'] not in bogus_values:
        attributes.append(['Package', biosample['package']])
    if biosample['submission_date'] not in bogus_values:
        attributes.append(['Submission date', biosample['submission_date']])
    if biosample['publication_date'] not in bogus_values:
        attributes.append(['Publication date', biosample['publication_date']])
    if biosample['last_updated'] not in bogus_values:
        attributes.append(['Last updated', biosample['last_updated']])

    # Sample Statistics
    assembly_info = result_detail.json()['reports'][0]['assembly_info']
    assembly_stats = result_detail.json()['reports'][0]['assembly_stats']

    if 'total_sequence_length' in assembly_stats:
        attributes.append(['Genome size', assembly_stats['total_sequence_length']])
    if 'total_number_of_chromosomes' in assembly_stats:
        attributes.append(['Number of chromosomes', assembly_stats['total_number_of_chromosomes']])
    if 'number_of_scaffolds' in assembly_stats:
        attributes.append(['Number of scaffolds', assembly_stats['number_of_scaffolds']])
    if 'scaffold_n50' in assembly_stats:
        attributes.append(['Scaffold N50', assembly_stats['scaffold_n50']])
    if 'scaffold_l50' in assembly_stats:
        attributes.append(['Scaffold L50', assembly_stats['scaffold_l50']])
    if 'number_of_contigs' in assembly_stats:
        attributes.append(['Number of contigs', assembly_stats['number_of_contigs']])
    if 'contig_n50' in assembly_stats:
        attributes.append(['Contig N50', assembly_stats['contig_n50']])
    if 'contig_l50' in assembly_stats:
        attributes.append(['Contig L50', assembly_stats['contig_l50']])
    if 'gc_percent' in assembly_stats:
        attributes.append(['GC percent', assembly_stats['gc_percent']])
    if 'assembly_level' in assembly_info:
        attributes.append(['Assembly level', assembly_info['assembly_level']])

    try:
        data.append(attributes)
    except:
        data_failed.append(accession)
        return


def sub_process(job, data, data_failed):
    for accession in job:
        subProcess = multiprocessing.Process(target=worker, args=(accession, data, data_failed))
        subProcess.start()


def collect_links(taxon_id, accession_list):
    request_payload = '{"filters":{"has_annotation":true,"assembly_source":"refseq","exclude_paired_reports":true,"assembly_version":"current"},"page_size":1000,"page_token":null,"returned_content":"COMPLETE","taxons":["' + str(
        taxon_id) + '"]}'
    result_page = requests.post("https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/dataset_report",
                                data=request_payload)
    for i in result_page.json()['reports']:
        accession_list.append(i['accession'])

    while 'next_page_token' in result_page.json():
        next_page_token = result_page.json()['next_page_token']
        request_payload = '{"filters":{"has_annotation":true,"assembly_source":"refseq","exclude_paired_reports":true,"assembly_version":"current"},"page_size":1000,"page_token":"' + next_page_token + '","returned_content":"COMPLETE","taxons":["' + str(
            taxon_id) + '"]}'
        result_page = requests.post("https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/dataset_report",
                                    data=request_payload)
        for i in result_page.json()['reports']:
            accession_list.append(i['accession'])
    return accession_list


def collect_data(job_pool, data, data_failed):
    for i, job in enumerate(job_pool):
        # print("Now working on batch", i, "：making", multi_proc_num, "requests in parallel...", end='')
        process = multiprocessing.Process(target=sub_process, args=(job, data, data_failed))
        process.start()
        process.join()


def watcher(total, list, list_2=[]):
    while True:
        time.sleep(0.2)
        num = len(list) + len(list_2)
        i = round((num / total) * 100)
        progress_bar(num, total, i)
        if num == total:
            print()
            break


def progress_bar(done, total, i):
    # print("\r", end="")
    print("\r" + "▋" * (i // 4) + " {}% ".format(i) + "(" + str(done) + "/" + str(total) + ")", end="",
          flush=True)
    # sys.stdout.flush()


def check_taxon(taxon_id):
    request_payload = '{"filters":{"has_annotation":true,"assembly_source":"refseq","exclude_paired_reports":true,"assembly_version":"current"},"page_size":1000,"page_token":null,"returned_content":"COMPLETE","taxons":["' + str(
        taxon_id) + '"]}'
    try:
        result_page = requests.post("https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/dataset_report",
                                    data=request_payload)
        return result_page.json()['total_count']
    except:
        return -1


if __name__ == '__main__':
    print('- NCBI Metadata Parser by Ryan -')
    # check taxon id
    taxon_valid = 0
    while not taxon_valid:
        try:
            taxon_id = int(input('\n' + '\x1b[93m' + 'Taxon ID: ' + '\x1b[0m'))
        except:
            print('\x1b[0;30;41m' + 'Syntax error.' + '\x1b[0m')
            continue
        print("Checking ID...")
        total_links_count = check_taxon(taxon_id)
        if total_links_count <= 0:
            print('\x1b[0;30;41m' + 'Bad ID.' + '\x1b[0m')
        else:
            print('\x1b[6;30;42m' + 'Legit ID.' + '\x1b[0m')
            taxon_valid = 1

    # multiprocessing config
    multi_config_valid = 0
    while not multi_config_valid:
        try:
            multi_proc_num = int(input('\n' + '\x1b[93m' + "Desired number of parallel processes: " + "\x1b[0m"))
            if multi_proc_num <= 0:
                print('\x1b[0;30;41m' + 'Bad number.' + '\x1b[0m')
                continue
            else:
                print('\x1b[6;30;42m' + 'Configuration saved.' + '\x1b[0m')
                multi_config_valid = 1
        except:
            print('\x1b[0;30;41m' + 'Syntax error.' + '\x1b[0m')
            continue

    manager = multiprocessing.Manager()
    accession_list = manager.list()

    # collect links
    print("\nFetching links:")
    process_collect_links = multiprocessing.Process(target=collect_links, args=(taxon_id, accession_list))
    process_watcher = multiprocessing.Process(target=watcher, args=(total_links_count, accession_list))
    process_collect_links.start()
    process_watcher.start()
    process_collect_links.join()
    process_watcher.join()
    print('\x1b[6;30;42m' + 'Done.' + '\x1b[0m')

    # allocate jobs
    print("\nAllocating jobs...")
    job_pool = []
    i = 0
    while i < len(accession_list):
        if (i + multi_proc_num - 1) <= (len(accession_list) - 1):
            job_pool.append(accession_list[i:(i + multi_proc_num)])
            i += multi_proc_num
        else:
            job_pool.append(accession_list[i:])
            i += multi_proc_num
    print('\x1b[6;30;42m' + 'Done.' + '\x1b[0m')

    # collect metadata
    print("\nCollecting and parsing data:")
    data = manager.list()
    data_failed = manager.list()
    process_watcher = multiprocessing.Process(target=watcher, args=(total_links_count, data, data_failed))
    process_collect_data = multiprocessing.Process(target=collect_data, args=(job_pool, data, data_failed))
    process_collect_data.start()
    process_watcher.start()
    process_collect_data.join()
    process_watcher.join()

    # retry
    while len(data_failed) > 0:
        retry = 1
        print("Number of failed requests:", len(data_failed))
        retry_valid = 0
        while not retry_valid:
            retry_confirmation = input('\x1b[93m' + "Retry them? (yes/no): " + "\x1b[0m")
            if retry_confirmation == 'yes' or retry_confirmation == 'no':
                if retry_confirmation == 'yes':
                    retry_valid = 1
                    print('\x1b[6;30;42m' + 'OK.' + '\x1b[0m')
                else:
                    retry_valid = 1
                    retry = 0
                    print('\x1b[6;30;42m' + 'OK.' + '\x1b[0m')
            else:
                retry_valid = 0
                print('\x1b[0;30;41m' + 'Syntax error.' + '\x1b[0m')
                continue

        if retry == 0:
            break

        print("\nRetrying failed links...")
        total_links_count_failed = len(data_failed)
        job_pool = []
        i = 0
        while i < len(data_failed):
            if (i + multi_proc_num - 1) <= (len(data_failed) - 1):
                job_pool.append(data_failed[i:(i + multi_proc_num)])
                i += multi_proc_num
            else:
                job_pool.append(data_failed[i:])
                i += multi_proc_num

        data_failed = manager.list()
        process_watcher = multiprocessing.Process(target=watcher, args=(total_links_count_failed, data_failed))
        process_collect_data = multiprocessing.Process(target=collect_data, args=(job_pool, data, data_failed))
        process_collect_data.start()
        process_watcher.start()
        process_collect_data.join()
        process_watcher.join()

    print('\x1b[6;30;42m' + 'Done.' + '\x1b[0m')

    # Sorting
    print('\nSorting data...')
    headers = []
    for item in data:
        for sub_item in item:
            if sub_item[0] not in headers:
                headers.append(sub_item[0])

    data_sorted = []
    for item in data:
        headers_raw = []
        data_raw = []
        data_sorted_single = []
        for sub_item in item:
            headers_raw.append(sub_item[0])
            data_raw.append(sub_item[1])
        for header in headers:
            if header in headers_raw:
                data_sorted_single.append(data_raw[headers_raw.index(header)])
            else:
                data_sorted_single.append('')
        data_sorted.append(data_sorted_single)
    print('\x1b[6;30;42m' + 'Done.' + '\x1b[0m')

    # output
    print('\nOutputting to csv...')
    with open('result_taxonid_' + str(taxon_id) + '_time_' + str(time.time()) + '.csv', 'w+') as out:
        writer = csv.writer(out)
        writer.writerow(headers)
        for row in data_sorted:
            writer.writerow(row)
    print('\x1b[6;30;42m' + 'Done.' + '\x1b[0m')
