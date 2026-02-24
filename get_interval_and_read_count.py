import pandas as pd

inputfile = "m6A_pos_on_genome_xgb_HEK293T_WT3.tsv"
outputfile = "intersections_m6A_sites_xgb_HEK293T_WT3.tsv"


def readintsv(filename, field_to_group):
    df = pd.read_csv(filename, sep = "\t", header=0)
    df_split = df.groupby(field_to_group)
    return df_split





def find_interval_intersections(intervals):
    # Sort intervals based on start positions
    sorted_intervals = sorted(intervals, key=lambda x: (x[0], x[1], x[2]))

    intersections = []
    current_intersection = sorted_intervals[0]
    unique_reads = set([current_intersection[2]])
    readname_count_temp_obj = {current_intersection[2]:1}
    count_intervals = 1
    for interval in sorted_intervals[1:]:
        # Check for intersection
        if current_intersection[1] > interval[0]:
            print("interval",interval)
            print("intersection",current_intersection)
            count_intervals += 1
            if interval[2] == current_intersection[2]:
                current_intersection = (max(current_intersection[0], interval[0]), min(current_intersection[1], interval[1]), interval[2])
                if interval[2] in readname_count_temp_obj:
                    readname_count_temp_obj[interval[2]] = readname_count_temp_obj[interval[2]] + 1
                else:
                    readname_count_temp_obj[interval[2]] = 1
            else:
                readname_count_temp_obj[interval[2]] = 1
                current_intersection = (max(current_intersection[0],interval[0]),min(current_intersection[1],interval[1]), interval[2])
            unique_reads.add(interval[2])

        else:
            # Save the current intersection and reset
            readname_count_array = [readname_count_temp_obj[key] for key in readname_count_temp_obj.keys()]

            intersections.append((current_intersection, count_intervals, len(unique_reads), readname_count_array))
            current_intersection = interval
            readname_count_temp_obj = {current_intersection[2]:1}
            count_intervals = 1
            unique_reads = set([current_intersection[2]])

    # Save the last intersection
    readname_count_array = [readname_count_temp_obj[key] for key in readname_count_temp_obj.keys()]
    intersections.append((current_intersection, count_intervals, len(unique_reads), readname_count_array))

    return intersections

def main():
    df_split = readintsv(inputfile,'Chromosome')
    intersection = {}
    for key,df_value in df_split:
        chrom_cur = key
        df_value_by_flag = df_value.groupby("Flag")
        intersection_cur = {}
        for key_flag, df_flag in df_value_by_flag:
            interval = []

            for index in df_flag.index.values:
                interval_cur = (df_flag["Start"][index],df_flag["End"][index],df_flag['Readname'][index])
                interval.append(interval_cur)
            intersection_cur[key_flag] = find_interval_intersections(interval)
        intersection[chrom_cur] = intersection_cur

    chrom = []
    startpos = []
    endpos = []
    count_in = []
    count_r = []
    flag = []
    readname_count_array = []
    for i in intersection.keys():
        for j in intersection[i].keys():
            for interval, count_interval, count_reads, readname_count in intersection[i][j]:
                chrom.append(i)
                startpos.append(interval[0])
                endpos.append(interval[1])
                count_in.append(count_interval)
                count_r.append(count_reads)
                flag.append(j)
                readname_count_array.append(readname_count)
    intersection_df = pd.DataFrame({"Chromosome":chrom, "Start":startpos, "End": endpos, "count_interval": count_in, "count_reads":count_r, "Flag": flag, "Count_Read_Interval": readname_count_array})
    intersection_df.to_csv(outputfile,sep='\t')
    return

if __name__ == '__main__':
    main()
