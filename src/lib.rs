use std::cmp::Reverse;
use std::collections::{BinaryHeap, BTreeMap};
use std::iter::FromIterator;
use hashbrown::{HashSet, HashMap};
use log::{info};
use smallvec::SmallVec;


//--------------------------------OLD--------------------------------------------










/// Interval vector sorting
pub fn sort_vector(intervals: &mut Vec<(u32, u32)>){
    intervals.sort_by(|a, b| (a.0.cmp(&b.0).then(b.1.cmp(&a.1))));

}



//--------------------------------------NEW---------------------------------------------------------------


/// Create which intervals are overlapping
/// Input:
/// - intervals_sorted: & Vec<(u32, u32)>
/// Output:
/// - Vec(Tuple to tuple)
pub fn make_nested2<'a, 'b>(intervals_sorted: &'a Vec<[u32; 2]>) -> Vec<(&[u32; 2], &[u32; 2])>{
    //info!("dasjkdjada");
    //let mut itervals = Vec::with_capacity(intervals_sorted.len());
    let mut open_intervals: BinaryHeap<Reverse<&[u32; 2]>> = BinaryHeap::new();

    // The result vector
    let mut result = Vec::with_capacity(intervals_sorted.len());
    // Iterate over sorted unique interval vector
    for tuple1 in intervals_sorted.iter() {
        // Remove tuples which have a smaller end then the start of the new tuple
        open_intervals.retain(|a| tuple1[0] < a.0[1]);

        // If empty (just at the start), add to open-intervals
        if !open_intervals.is_empty() {

            // Sort the open intervals by start and end


            // Get those that have bigger end than the one you are looking at
            for interval in open_intervals.iter(){
                if interval.0[1] >= tuple1[1] - 1{
                    result.push((tuple1, interval.0));
                } else {
                    break
                }
            }

            // println!("open {:?}", open_intervals);
            // println!("tuple {:?}", tuple1);
            // println!("result {:?}\n", result);
        }
        open_intervals.push(Reverse(tuple1));

    }
    result.shrink_to_fit();
    return result
}


/// Create which intervals are overlapping
/// Input:
/// - intervals_sorted: & Vec<(u32, u32)>
/// Output:
/// - Vec(Tuple to tuple)
pub fn make_nested22<'a, 'b>(intervals_sorted: &'a Vec<(u32, u32)>) -> HashMap<&(u32, u32), Vec<&(u32, u32)>>{
    //info!("dasjkdjada");
    //let mut itervals = Vec::with_capacity(intervals_sorted.len());
    let mut open_intervals: BinaryHeap<Reverse<&(u32, u32)>> = BinaryHeap::new();

    // The result vector
    let mut result = HashMap::with_capacity(intervals_sorted.len());
    // Iterate over sorted unique interval vector
    for tuple1 in intervals_sorted.iter() {
        // Remove tuples which have a smaller end then the start of the new tuple
        open_intervals.retain(|a| tuple1.0 < a.0.1);

        // If empty (just at the start), add to open-intervals
        if !open_intervals.is_empty() {

            // Sort the open intervals by start and end


            //Get those that have bigger end than the one you are looking at
            for interval in open_intervals.iter(){
                if interval.0.1 >= tuple1.1 - 1{
                    result.entry(tuple1).or_insert(vec![]).push( interval.0)
                } else {
                    break
                }
            }

            // println!("open {:?}", open_intervals);
            // println!("tuple {:?}", tuple1);
            // println!("result {:?}\n", result);
        }
        open_intervals.push(Reverse(tuple1));

    }
    result.shrink_to_fit();
    return result
}


//-----------------Tests----------------------------------------------------------------------


#[cfg(test)]
mod tests {
    use log::info;
    use crate::{sort_vector, make_nested2};


    fn init() {
        env_logger::Builder::new().filter_level(log::LevelFilter::Trace).try_init().err();
    }
    // cargo test -- --nocapture




    #[test]
    fn basic() {
        init();
        // We test remove and and general function
        info!("\nRunning test 1");
        //assert_eq!(2 + 2, 4);

        let mut intervals: Vec<(u32, u32)> = vec![(1,10), (2,4), (7,9), (8,9), (1,20), (3,20), (1,3)];;
        assert_eq!(intervals.len(), 7);
        sort_vector(&mut intervals);
        info!("");
    }

}
