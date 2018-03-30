#!/usr/bin/python3
import pandas
import os

def main():
    sheet_name = 'Slope_match_PacBio_cells_Diane_March23_2018.xlsx'
    clustergram_sheet = pandas.read_excel(sheet_name, sheetname='clustergram', index_col='cell_name')
    slope_sheet = pandas.read_excel(sheet_name, sheetname='slope',
                              index_col='cell_name'
    )

    slope = {}
    for i, row in slope_sheet.iterrows():
        key = row.name.replace('_mm10','').replace('_clean', '')
        if not key.endswith('.1'):
            slope[key] = row.values[0]

    slope_col = []
    for i, row in clustergram_sheet.iterrows():
        slope_col.append(slope.get(row.name))

    clustergram_sheet['slope'] = slope_col
    print(clustergram_sheet)

    missing_from_clustergram = set(clustergram_sheet.index).difference(set(slope.keys()))
    print('missing clustergram', sorted(missing_from_clustergram))
    missing_from_slope = set(slope.keys()).difference(set(clustergram_sheet.index))
    print('missing slope', sorted(missing_from_slope))
    clustergram_sheet.to_csv('clustergram_slope_merged.csv')
    
if __name__ == "__main__":
    main()
