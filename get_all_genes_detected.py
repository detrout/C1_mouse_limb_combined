#!/usr/bin/python3

import psycopg2

def main():
    out = get_bams_by_library()

def get_bams_by_library():
    conn = psycopg2.connect("dbname=htsworkflow")
    cur = conn.cursor()

    cur.execute("""
with
  experiment as (
  select uri as Experiment,
         payload->>'accession' as Experiment_Accession,
         payload->>'description' as Experiment_Description,
         payload->>'status' as Experiment_Status,
         payload->>'date_released' as Experiment_Released,
         payload->>'assay_title' as Experiment_Type,
         jsonb_array_elements_text(payload->'replicates') as Replicate
  from item
  where object_type = 'Experiment'
  ),
  replicate as (
    select uri as Replicate,
           payload->>'library' as Library,
           payload->>'antibody' as Antibody
    from item
    where object_type = 'Replicate'
  ),
  library as (
    select uri as Library,
           payload->>'accession' as Library_Accession,
           payload->>'date_created' as Library_Created,
           payload->>'biosample' as Biosample
    from item
    where object_type = 'Library'
  )
select Experiment, library.Library
from experiment
  left join replicate on experiment.replicate = replicate.Replicate
  left join library on replicate.Library = library.Library
limit 10
""")
    print(cur.fetchone())

if __name__ == '__main__':
    main()
