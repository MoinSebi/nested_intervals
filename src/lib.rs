use std::cmp::Reverse;
use std::collections::{BinaryHeap, BTreeMap};
use std::iter::FromIterator;
use hashbrown::{HashSet, HashMap};
use log::{info};
use smallvec::SmallVec;


/// Create which intervals are overlapping
/// Input:
/// - intervals_sorted: & Vec<(u32, u32)>
/// Output:
/// - Vec(Tuple to tuple)
pub fn make_nested<'a, 'b>(intervals_sorted: &'a Vec<[u32; 3]>) -> Vec<(u32, u32)>{
    //info!("dasjkdjada");
    //let mut itervals = Vec::with_capacity(intervals_sorted.len());
    let mut open_intervals: Vec<&[u32; 3]> = Vec::new();
    let mut result = Vec::with_capacity(intervals_sorted.len());
    // Iterate over sorted unique interval vector
    let mut d = Vec::new();
    for tuple1 in intervals_sorted.iter() {
        open_intervals.retain(|a| tuple1[0] < a[1]);

        // If empty (just at the start), add to openintervals
        if open_intervals.len() != 0 {
            open_intervals.sort_by_key(|a| (Reverse(a[0]), a[1]));



            if let Some(index) = open_intervals.iter().position(|a| a[1] >= tuple1[1]){
                let mut smallest = &open_intervals.get(index).unwrap()[1];
                result.push((tuple1[2], open_intervals.get(index).unwrap()[2]));
                d = open_intervals.iter().filter(|a| smallest >  &a[1]).collect();
                //println!("{:?}", d);
                // println!("{:?}", open_intervals);
                // println!("{:?}", smallest);
                // println!("{:?}\n", tuple1);
                for x in open_intervals.iter() {
                    if &x[1] < smallest && x[1] > tuple1[1] {
                        smallest = &x[1];
                        result.push((tuple1[2], x[2]));
                    }
                }
            }
        }
        open_intervals.push(tuple1);
    }
    result.shrink_to_fit();
    return result
}

//-----------------Helper----------------------------------------------------------------------

/// Sort the interval by first entry and reverse by second
///
/// ```
///
/// ```
pub fn sort_vector(intervals: &mut Vec<[u32; 3]>){
    intervals.sort_by(|a, b| (a[0].cmp(&b[0]).then(b[1].cmp(&a[1]))));

}


//-----------------Tests----------------------------------------------------------------------


#[cfg(test)]
mod tests {
    use log::info;
    use crate::{make_nested, sort_vector};


    fn init() {
        env_logger::Builder::new().filter_level(log::LevelFilter::Trace).try_init().err();
    }
    // cargo test -- --nocapture




    #[test]
    fn basic()  {
        init();
        // We test remove and and general function
        info!("\nRunning test 1");
        //assert_eq!(2 + 2, 4);

        let mut intervals: Vec<[u32; 3]> = vec![[1,10,1], [2,4,2], [7,9,3], [8,9,4], [1,20,5], [3,20,6], [1,3, 7]];;
        assert_eq!(intervals.len(), 7);
        sort_vector(&mut intervals);
        let intersection = make_nested(&intervals);
        let vecc: Vec<(u32, u32)> = Vec::new();
        assert_eq!(intersection, vecc);
    }

}
