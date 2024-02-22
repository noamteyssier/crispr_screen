use crate::aggregation::GeneAggregation;
use anyhow::{bail, Result};
use csv::ReaderBuilder;
use ndarray::{Array1, Array2, Axis};
use std::collections::HashMap;
use std::{
    fs::File,
    io::{BufRead, BufReader},
};

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
        Ok(File::open(filename).map(BufReader::new)?)
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
            map.entry(s.clone()).or_default();
            map
        })
    }

    fn build_data_hashmap(headers: &[String]) -> HashMap<String, Vec<f64>> {
        headers.iter().skip(2).fold(HashMap::new(), |mut map, s| {
            map.entry(s.clone()).or_default();
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
                    let value = record
                        .get(h)
                        .unwrap_or_else(|| panic!("Malformed Record, missing header: {h}"));
                    meta_hash.get_mut(h).unwrap().push(value.to_string());
                } else {
                    let value = match record.get(h) {
                        Some(x) => x.parse::<f64>().unwrap_or(0.0),
                        None => bail!("Malformed Record, missing header: {h}"),
                    };
                    data_hash.get_mut(h).unwrap().push(value);
                }
            }
        }

        Self::validate_input(headers, meta_hash, data_hash)
    }

    fn validate_input(
        headers: &[String],
        meta_hash: &HashMap<String, Vec<String>>,
        data_hash: &HashMap<String, Vec<f64>>,
    ) -> Result<()> {
        // Checks if there were not enough headers read
        if headers.len() != meta_hash.len() + data_hash.len() {
            bail!("Malformed data entry in Axis(0)")
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
        if labels.len() > self.data.keys().len() {
            bail!("Provided count matrix is of unexpected shape!\nMust be [sgrna, gene, samples ... ]")
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

    pub fn get_f64_column(&self, name: &str) -> Result<Array1<f64>> {
        if let Some(column) = self.data.get(name) {
            Ok(Array1::from_vec(column.clone()))
        } else {
            bail!("Column not found")
        }
    }

    pub fn get_string_column(&self, name: &str) -> Result<Vec<String>> {
        if let Some(column) = self.meta.get(name) {
            Ok(column.clone())
        } else {
            bail!("Column not found")
        }
    }

    pub fn validate_ntc(&self, config: &GeneAggregation) -> Result<()> {
        match config {
            GeneAggregation::Inc {
                token,
                fdr: _,
                group_size: _,
                n_draws: _,
                use_product: _,
            } => {
                if self.get_sgrna_names().iter().any(|x| x.contains(token)) {
                    Ok(())
                } else {
                    bail!("Non-Targeting Token ({token}) not found in any sgrna names - please use RRA or update the provided token.")
                }
            }
            _ => Ok(()),
        }
    }
}

#[cfg(test)]
mod testing {
    use crate::aggregation::GeneAggregation;

    use super::SimpleFrame;
    use ndarray::s;
    use ndarray_rand::rand::random;

    fn example_dataset() -> String {
        let mut s = String::with_capacity(1000);
        let headers = ["sgrna", "gene", "low_1", "low_2", "high_1", "high_2"];
        let num_rows = 10;

        headers.iter().enumerate().for_each(|(idx, x)| {
            if idx > 0 {
                s.push_str(&format!("\t{x}"))
            } else {
                s.push_str(x)
            }
        });
        s.push('\n');

        (0..num_rows).for_each(|row_id| {
            headers.iter().enumerate().for_each(|(idx, _)| {
                if idx == 0 {
                    s.push_str(&format!("sgrna_{row_id}"));
                } else if idx == 1 {
                    s.push_str(&format!("\tgene_{}", row_id % 5));
                } else {
                    s.push_str(&format!("\t{}", random::<f64>()))
                }
            });
            s.push('\n');
        });
        s
    }

    fn broken_dataset() -> String {
        let mut s = String::with_capacity(1000);
        let headers = ["sgrna", "low_1", "low_2", "high_1", "high_2"];
        let num_rows = 10;

        headers.iter().enumerate().for_each(|(idx, x)| {
            if idx > 0 {
                s.push_str(&format!("\t{x}"))
            } else {
                s.push_str(x)
            }
        });
        s.push('\n');

        (0..num_rows).for_each(|row_id| {
            headers.iter().enumerate().for_each(|(idx, _)| {
                if idx == 0 {
                    s.push_str(&format!("sgrna_{row_id}"));
                } else {
                    s.push_str(&format!("\t{}", random::<f64>()))
                }
            });
            s.push('\n');
        });
        s
    }

    fn malformed_cells_empty() -> String {
        let mut s = String::with_capacity(1000);
        let headers = ["sgrna", "gene", "low", "high"];
        let num_rows = 10;

        // add headers
        headers.iter().enumerate().for_each(|(idx, x)| {
            if idx > 0 {
                s.push_str(&format!("\t{x}"))
            } else {
                s.push_str(x)
            }
        });
        s.push('\n');

        // add data
        (0..num_rows).for_each(|row_id| {
            headers.iter().enumerate().for_each(|(idx, _)| {
                if idx == 0 {
                    s.push_str(&format!("sgrna_{row_id}"));
                } else {
                    let n = random::<f64>();
                    if n > 0.1 {
                        s.push_str(&format!("\t{}", n))
                    } else {
                        s.push('\t')
                    }
                }
            });
            s.push('\n');
        });

        s
    }

    fn malformed_cells_string() -> String {
        let mut s = String::with_capacity(1000);
        let headers = ["sgrna", "gene", "low", "high"];
        let num_rows = 10;

        // add headers
        headers.iter().enumerate().for_each(|(idx, x)| {
            if idx > 0 {
                s.push_str(&format!("\t{x}"))
            } else {
                s.push_str(x)
            }
        });
        s.push('\n');

        // add data
        (0..num_rows).for_each(|row_id| {
            headers.iter().enumerate().for_each(|(idx, _)| {
                if idx == 0 {
                    s.push_str(&format!("sgrna_{row_id}"));
                } else {
                    let n = random::<f64>();
                    if n > 0.1 {
                        s.push_str(&format!("\t{}", n))
                    } else {
                        s.push_str("\tNaN")
                    }
                }
            });
            s.push('\n');
        });

        s
    }

    #[test]
    fn test_simple_frame() {
        let datastream = example_dataset();
        let frame = SimpleFrame::from_string(&datastream).unwrap();

        assert_eq!(
            frame.headers,
            vec!["sgrna", "gene", "low_1", "low_2", "high_1", "high_2"]
        );

        let dm = frame
            .data_matrix(
                &vec!["low_1", "low_2", "high_1", "high_2"]
                    .into_iter()
                    .map(|x| x.to_owned())
                    .collect::<Vec<String>>(),
            )
            .unwrap();
        assert_eq!(dm.shape(), &[10, 4]);
    }

    #[test]
    fn test_simple_frame_subset() {
        let datastream = example_dataset();
        let frame = SimpleFrame::from_string(&datastream).unwrap();

        let dm_sub = frame
            .data_matrix(
                &vec!["low_1", "high_1"]
                    .into_iter()
                    .map(|x| x.to_owned())
                    .collect::<Vec<String>>(),
            )
            .unwrap();
        assert_eq!(dm_sub.shape(), &[10, 2]);
    }

    #[test]
    fn test_simple_frame_ordering() {
        let datastream = example_dataset();
        let frame = SimpleFrame::from_string(&datastream).unwrap();

        let dm_rev = frame
            .data_matrix(
                &vec!["high_1", "low_1"]
                    .into_iter()
                    .map(|x| x.to_owned())
                    .collect::<Vec<String>>(),
            )
            .unwrap();
        let dm_fwd = frame
            .data_matrix(
                &vec!["low_1", "high_1"]
                    .into_iter()
                    .map(|x| x.to_owned())
                    .collect::<Vec<String>>(),
            )
            .unwrap();
        assert_eq!(dm_rev.slice(s![.., 0]), dm_fwd.slice(s![.., 1]));
        assert_eq!(dm_rev.slice(s![.., 1]), dm_fwd.slice(s![.., 0]));
    }

    #[test]
    fn test_missing_label() {
        let datastream = example_dataset();
        let frame = SimpleFrame::from_string(&datastream).unwrap();
        let dm = frame.data_matrix(
            &vec!["low_1", "low_2", "high_1", "missing_label"]
                .into_iter()
                .map(|x| x.to_owned())
                .collect::<Vec<String>>(),
        );
        assert!(dm.is_err());
    }

    #[test]
    fn test_broken_dataset() {
        let datastream = broken_dataset();
        let frame = SimpleFrame::from_string(&datastream).unwrap();
        assert_eq!(
            frame.headers,
            vec!["sgrna", "low_1", "low_2", "high_1", "high_2"]
        );
        let dm = frame.data_matrix(
            &vec!["low_1", "low_2", "high_1", "high_2"]
                .into_iter()
                .map(|x| x.to_owned())
                .collect::<Vec<String>>(),
        );
        assert!(dm.is_err());
    }

    #[test]
    fn test_gene_names() {
        let datastream = example_dataset();
        let frame = SimpleFrame::from_string(&datastream).unwrap();
        assert_eq!(
            frame.get_gene_names(),
            &vec![
                "gene_0", "gene_1", "gene_2", "gene_3", "gene_4", "gene_0", "gene_1", "gene_2",
                "gene_3", "gene_4"
            ]
        );
    }

    #[test]
    fn test_sgrna_names() {
        let datastream = example_dataset();
        let frame = SimpleFrame::from_string(&datastream).unwrap();
        assert_eq!(
            frame.get_sgrna_names(),
            &vec![
                "sgrna_0", "sgrna_1", "sgrna_2", "sgrna_3", "sgrna_4", "sgrna_5", "sgrna_6",
                "sgrna_7", "sgrna_8", "sgrna_9"
            ]
        );
    }

    #[test]
    fn test_malformed_cells_empty() {
        let datastream = malformed_cells_empty();
        let frame = SimpleFrame::from_string(&datastream).unwrap();
        assert_eq!(frame.headers, vec!["sgrna", "gene", "low", "high"]);
        let dm = frame.data_matrix(
            &vec!["low", "high"]
                .into_iter()
                .map(|x| x.to_owned())
                .collect::<Vec<String>>(),
        );
        assert!(dm.is_ok());
    }

    #[test]
    fn test_malformed_cells_string() {
        let datastream = malformed_cells_string();
        let frame = SimpleFrame::from_string(&datastream).unwrap();
        assert_eq!(frame.headers, vec!["sgrna", "gene", "low", "high"]);
        let dm = frame.data_matrix(
            &vec!["low", "high"]
                .into_iter()
                .map(|x| x.to_owned())
                .collect::<Vec<String>>(),
        );
        assert!(dm.is_ok());
    }

    #[test]
    fn test_ntc_token_validation() {
        let datastream = example_dataset();
        let frame = SimpleFrame::from_string(&datastream).unwrap();

        // Expected to be missing
        let config = GeneAggregation::Inc {
            token: "non-targeting",
            fdr: 0.05,
            group_size: 3,
            n_draws: 100,
            use_product: false,
        };
        let result = frame.validate_ntc(&config);
        assert!(result.is_err());

        // Expected to be found
        let config = GeneAggregation::Inc {
            token: "sgrna_2",
            fdr: 0.05,
            group_size: 3,
            n_draws: 100,
            use_product: false,
        };
        let result = frame.validate_ntc(&config);
        assert!(result.is_ok());
    }
}
