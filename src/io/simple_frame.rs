use std::{
    fs::File,
    io::{BufRead, BufReader},
};

use anyhow::{bail, Result};
use csv::ReaderBuilder;
use ndarray::{Array1, Array2, Axis};
use std::collections::HashMap;

type Record = HashMap<String, String>;

pub struct SimpleFrame {
    headers: Vec<String>,
    meta: HashMap<String, Vec<String>>,
    data: HashMap<String, Vec<f64>>,
}
impl SimpleFrame {
    pub fn from_filepath(path: &str) -> Result<Self> {
        let headers = Self::parse_headers(&mut Self::read_file(path)?)?;
        let mut meta = Self::build_meta_hashmap(&headers);
        let mut data = Self::build_data_hashmap(&headers);
        Self::parse_to_hashmaps(&mut Self::read_file(path)?, &headers, &mut meta, &mut data)?;
        Ok(Self {
            headers,
            meta,
            data,
        })
    }

    #[allow(dead_code)]
    pub fn from_string(string: &str) -> Result<Self> {
        let headers = Self::parse_headers(&mut string.as_bytes())?;
        let mut meta = Self::build_meta_hashmap(&headers);
        let mut data = Self::build_data_hashmap(&headers);
        Self::parse_to_hashmaps(&mut string.as_bytes(), &headers, &mut meta, &mut data)?;
        Ok(Self {
            headers,
            meta,
            data,
        })
    }

    fn read_file(filename: &str) -> Result<BufReader<File>> {
        Ok(
            File::open(filename)
                .map(BufReader::new)?
        )
    }

    fn parse_headers<R: BufRead>(buffer: &mut R) -> Result<Vec<String>> {
        let mut header_row = String::new();
        buffer.read_line(&mut header_row)?;
        let headers = header_row
            .trim()
            .split('\t')
            .map(std::string::ToString::to_string)
            .collect::<Vec<String>>();
        Ok(headers)
    }

    fn build_meta_hashmap(headers: &[String]) -> HashMap<String, Vec<String>> {
        headers.iter().take(2).fold(HashMap::new(), |mut map, s| {
            map.entry(s.clone()).or_insert(Vec::new());
            map
        })
    }

    fn build_data_hashmap(headers: &[String]) -> HashMap<String, Vec<f64>> {
        headers.iter().skip(2).fold(HashMap::new(), |mut map, s| {
            map.entry(s.clone()).or_insert(Vec::new());
            map
        })
    }

    fn parse_to_hashmaps<R: BufRead>(
        buffer: &mut R,
        headers: &[String],
        meta_hash: &mut HashMap<String, Vec<String>>,
        data_hash: &mut HashMap<String, Vec<f64>>,
    ) -> Result<()> {
        let mut reader = ReaderBuilder::new()
            .has_headers(true)
            .delimiter(b'\t')
            .from_reader(buffer);

        for result in reader.deserialize() {
            let record: Record = result?;
            for (idx, h) in headers.iter().enumerate() {
                if idx < 2 {
                    meta_hash
                        .get_mut(h)
                        .unwrap()
                        .push(record.get(h).expect("Malformed Record").clone());
                } else {
                    data_hash.get_mut(h).unwrap().push(
                        record
                            .get(h)
                            .expect("Malformed Record")
                            .parse::<f64>()
                            .expect("Unable to parse record to float")
                            .to_owned(),
                    );
                }
            }
        }
        Ok(())
    }

    pub fn valid_headers(&self, labels: &[String]) -> bool {
        labels.iter().all(|x| self.headers.contains(x))
    }

    pub fn data_matrix(&self, labels: &[String]) -> Result<Array2<f64>> {
        if !self.valid_headers(labels) {
            bail!("Invalid headers provided")
        }
        let matrix = labels
            .iter()
            .map(|x| {
                Array1::from_iter(
                    self.data
                        .get(x)
                        .expect("Unexpected column name provided")
                        .iter()
                        .copied(),
                )
                .insert_axis(Axis(0))
            })
            .reduce(|mut x, y| {
                x.push(Axis(0), y.view().remove_axis(Axis(0)))
                    .expect("Could not append column");
                x
            })
            .expect("Could not create matrix");
        let matrix = matrix.t().to_owned();
        Ok(matrix)
    }

    pub fn get_sgrna_names(&self) -> &[String] {
        self.meta.get(&self.headers[0]).unwrap()
    }

    pub fn get_gene_names(&self) -> &[String] {
        self.meta.get(&self.headers[1]).unwrap()
    }
}
