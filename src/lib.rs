use std::cmp::Reverse;
use std::collections::BTreeMap;
use std::iter::FromIterator;
use hashbrown::{HashSet, HashMap};
use log::{info};
use smallvec::SmallVec;


//--------------------------------OLD--------------------------------------------
/// Gives parent and child relation ship
/// Make vector and no hashset
/// get some more values or just delete is at all
/// use node directly 
#[derive(Eq, PartialEq, Debug, Clone)]
pub struct Network<'a>{
    pub parent: HashSet<&'a (u32, u32)>,
}


impl Network<'_> {
    /// Create empty network
    pub fn new() -> Self {
        let parent: HashSet<&(u32, u32)> = HashSet::new();
        Self {
            parent: parent,
        }
    }
}


/// Fill the graph structure
/// Input:
///     - intervals_sorted Vec<(u32, u32)>: Sorted interval (no duplicates)
///     - order: HashMap<(u32, u32), Network>): Interval relation (this is empty)
///
/// Output:
///     -
///
/// Running:
///     - Iterate over all intervals and sort them to open bubbles
///     - If there is some kind of relation -> Add to parent, add to child
///     - Some kind of wrapper function
/// Check:
///     - checker_rev
pub fn make_nested<'a, 'b>(intervals_sorted: &'a Vec<(u32, u32)>, network: &'b mut BTreeMap<&'a (u32, u32), Network<'a>>){
    let mut open_intervals = SmallVec::<[&(u32, u32); 20]>::new();
    
    // Iterate over sorted unique interval vector
    for tuple1 in intervals_sorted.iter(){
        // If empty (just at the start), add to openintervals
        if open_intervals.len() != 0{

            let mut hits = SmallVec::<[&(u32, u32); 20]>::new();
            let mut overlaps = SmallVec::<[&(u32, u32); 20]>::new();

            // Iterate over open intervals
            // For each interval: Start checking, when and if hit or overlap
            for open1 in open_intervals.iter(){
                let mut hs = Vec::new();
                let solution = checker_recursive(open1, tuple1, network, false, & mut hs);
                hits.extend(solution.0);
                overlaps.extend(solution.1);
            }
            if hits.len() == 0{
                //info!("Nothing happening");
            }
            else {
                //info!("We need to filter");


                filter_hit(& mut hits);
                let hs: HashSet<&(u32, u32)> = hits.into_iter().collect();
                for x in hs{
                    network.get_mut(tuple1).unwrap().parent.insert(x);
                }


            }

            open_intervals =  SmallVec::<[&(u32, u32); 20]>::new();
            for x in overlaps.into_iter() {
                open_intervals.push(x);
            }
            open_intervals.push(tuple1);

        } else {
            open_intervals.push(&tuple1)


        }
       //info!("OP {:?}", &open_intervals);
        //info!("{:?}", order);

    }

    //info!("{:?}", order);
}



/// In a list of intervals, check if some are related to each other (very slow)
/// Related: -> One is parent/child of another
/// Complexity: O(n^2)
/// Input:
///     - candidates Vec<(u32, u32)>: List of intervals to check
///
/// Mut:
///     - candidates Vec<(u32, u32)>: List of intervals to check
///
/// Output:
///     --> Chaneged candidates
///
/// Running:
///     - Iterate over the list, check if the relation between the candidates
///
/// TODO:
/// - Stuff is removed at the end, checking like everything
/// PROBLEM:
/// - O(n^2) complexity
pub fn filter_hit(candidates: &mut SmallVec<[&(u32, u32); 20]>) {
    let cand_vector: Vec<&(u32, u32)> = candidates.iter().cloned().collect();
    for (i1, x) in cand_vector.iter().enumerate(){
        for y in cand_vector[i1+1..].iter(){
            // x is außerhalb von y (oder andersrum)
            // THIS IS NOT TRUE LOL
            if (x.0 <= y.0) & (x.1 >= y.1) {
                let index = candidates.iter().position(|x1| x1 == x).unwrap();
                candidates.remove(index);
                //break;
                // x is der parent
            } else if (x.0 >= y.0) & (x.1 <= y.1){
                let index = candidates.iter().position(|x1| x1 == y).unwrap();
                candidates.remove(index);
            }
        }
    }
}







/// Recursive main function
/// Input:
///     - old: "old" interval
///     - new: "new" interval
///     - hm: main structure to change (needed for parents)
///     - overlaps_parent: bool
///     - hs: HashSet of all already tested intervals
///
/// Output:
///     - (hit, overlap)
///
/// Running:
///     1. Check if new interval is within the old one -> Hit and not overlap
///     2. if not (no hit or interval):
///         - if your parent is not overlapping and you overlap -> not interesting
///         - rerun everything with the parent (recursive)
///         - merge return
///
/// Problem:
///     - You can check the same interval multiple times
pub fn checker_recursive<'a, 'b>(old: &'a(u32, u32), new: &'a (u32, u32), hm: &'b mut BTreeMap<&'a ( u32, u32), Network<'a>>, overlaps_parent: bool, hs: & mut Vec<&'a (u32, u32)>) -> (Vec<&'a (u32, u32)>, Vec<&'a( u32, u32)>) {
    let mut hits: Vec<&(u32, u32)> = Vec::new();
    let mut overlaps: Vec<&(u32, u32)> = Vec::new();

    let mut now_overlapping: bool = false;

    // It is a hit!
    if (old.0 <= new.0) & (old.1 >= new.1) {
        hits.push(old);

        // not a hit
    } else {
        if (!overlaps_parent) & (old.1 > new.0) {
            overlaps.push(old);
            now_overlapping = true;
        }
        if hm.get(old).unwrap().parent.len() != 0 {
            let mut vecc_p = Vec::new();
            for x in hm.get(old).unwrap().parent.iter() {
                //u.append(&mut helper1(x, new, hm));
                vecc_p.push(*x);
            }
            for x in vecc_p.iter() {
                if !hs.contains(x){
                    hs.push(*x);
                    //info!("{:?}", x);
                    let jo = &mut checker_recursive(&x, new, hm, now_overlapping, hs);
                    //info!("{:?}", jo);
                    hits.extend(& jo.0);
                    overlaps.extend(& jo.1);
                }
            }
        }
    }
    (hits, overlaps)
}



/// Interval vector sorting
pub fn sort_vector(intervals: &mut Vec<(u32, u32)>){
    intervals.sort_by(|a, b| (a.0.cmp(&b.0).then(b.1.cmp(&a.1))));

}


/// Create the "network" for parent-child relationship
pub fn create_network_hashmap(intervals_sorted: & Vec<(u32, u32)>) -> BTreeMap<&(u32, u32), Network>{
    let mut order = BTreeMap::new();
    for tuple1 in intervals_sorted.iter(){
        order.insert(tuple1, Network::new());
    }
    order
}




//--------------------------------------NEW---------------------------------------------------------------

pub fn make_nested2<'a, 'b>(intervals_sorted: &'a Vec<(u32, u32)>) -> Vec<((u32, u32), &(u32, u32))>{
    //info!("dasjkdjada");
    //let mut itervals = Vec::with_capacity(intervals_sorted.len());
    let mut open_intervals: Vec<&(u32, u32)> = Vec::new();
    let mut result = Vec::with_capacity(intervals_sorted.len());
    // Iterate over sorted unique interval vector
    for tuple1 in intervals_sorted.iter() {
        open_intervals.retain(|a| tuple1.0 < a.1);

        // If empty (just at the start), add to openintervals
        if open_intervals.len() != 0 {
            open_intervals.sort_by_key(|a| (Reverse(a.0), a.1));


            let mut index = open_intervals.iter().position(|a| a.1 > tuple1.1).unwrap();
            //let mut index = open_intervals.partition_point(|a| a.1 > tuple1.1)-1;

            let mut smallest = &open_intervals.get(index).unwrap().1;
            result.push((*tuple1, *open_intervals.get(index).unwrap()));
            for x in open_intervals.iter() {
                if &x.1 < smallest && x.1 > tuple1.1 {
                    smallest = &x.1;
                    result.push((*tuple1, *x));
                }
            }
        }
        open_intervals.push(tuple1);
    }
    result.shrink_to_fit();
    return result
}



//-----------------HELPER----------------------------------------------------------------------

/// This is redundant to remove_duplicates
pub fn check_unique(intervals: & Vec<(u32, u32)>) -> bool{
    let hs: HashSet<(u32, u32)> = HashSet::from_iter(intervals.iter().cloned());
    if intervals.len() == hs.len(){
        true
    } else {
        false
    }
}


/// Removing duplicated intervals
///
/// If not removed, they are parent + child of the same interval
pub fn remove_duplicates(intervals: & mut Vec<(u32, u32)> ){
    let k2: HashSet<(u32, u32)> = intervals.iter().cloned().collect();
    let mut uniques = HashSet::new();
    if ! (k2.len() == intervals.len()){
        info!("Removed duplicates");
        intervals.retain(|e| uniques.insert(*e));
    }

}


/// Check if start and stop in order
/// Input:
///     - intervals: & Vec<(u32, u32)>)
/// Return:
///     - bool:
///         - false if not in order
///         - true if in order
/// Problem:
///     - not 100% right
pub fn start_stop_check(intervals: & Vec<(u32, u32)>) -> bool{
    for x in intervals.iter(){
        if x.0 < x.1{
            return false
        }
    }
    true
}


/// If entries are not in order, change them
pub fn start_stop_order(intervals: &mut Vec<(u32, u32)>){
    let mut change_list: Vec<usize> = Vec::new();
    for  (i, x)  in intervals.iter().enumerate(){
        if x.0 > x.1{
            change_list.push(i)

        }
    }
    for x in change_list.iter(){
        intervals[*x] = (intervals[*x].1, intervals[*x].0)
    }
}


/// Check if intervals are overlapping
/// If not, you can use the the much faster machine
pub fn check_overlapping(intervals: &mut Vec<(u32, u32)>) -> bool{
    for (i1, x) in intervals.iter().enumerate(){
        for y in intervals[i1+1..].iter(){
            if ((y.0 < x.0) & (x.0 < y.1) & (x.1 > y.1)) | ((x.0 < y.0) & (y.0 < x.1) & (y.1 > x.1)){
                return true

            }
        }
    }
    false
}




#[cfg(test)]
mod tests {
    use log::info;
    use crate::{sort_vector, make_nested, create_network_hashmap, remove_duplicates, make_nested2};


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

        let mut intervals: Vec<(u32, u32)> = vec![(1,10), (2,4), (7,9), (7,9), (8,9), (1,20)];
        remove_duplicates(&mut intervals);
        assert_eq!(intervals.len(), 5);
        sort_vector(&mut intervals);
        let mut network = create_network_hashmap(& intervals);

        make_nested(&intervals, & mut network);
        info!("{:?}", network);
        info!("");


    }

    #[test]
    fn basic2() {
        init();
        // We test remove and and general function
        info!("\nRunning test 1");
        //assert_eq!(2 + 2, 4);

        let mut intervals: Vec<(u32, u32)> = vec![(1,10), (2,4), (7,9), (7,9), (8,9), (1,20)];;
        remove_duplicates(&mut intervals);
        assert_eq!(intervals.len(), 5);
        sort_vector(&mut intervals);
        let f = make_nested2(&mut intervals);

        info!("this is f {:?}", f);
        info!("");
    }

    // #[test]
    // /// Overlap
    // fn run1(){
    //     init();
    //     info!("Test Nr. 2 -- Removing shit");
    //     let mut intervals: Vec<(u32, u32)> = vec![(1,20), (5,30)];
    //     sort_vector(&mut intervals);
    //     let mut network = create_network_hashmap(& intervals);
    //     make_nested(&intervals, & mut network);
    //
    //     assert_eq!(network.len(), 2);
    //     assert_eq!(intervals.len(), 2);
    // }
    //
    //
    // #[test]
    // /// double hit
    // fn run2(){
    //     init();
    //     info!("Test Nr. 2 -- Removing shit");
    //     let mut intervals: Vec<(u32, u32)> = vec![(1,20), (10,30), (12,18)];
    //     sort_vector(&mut intervals);
    //     let mut network = create_network_hashmap(& intervals);
    //     make_nested(&intervals, & mut network);
    //
    //     info!("run2{:?}", network);
    //     assert_eq!(intervals.len(), 3);
    // }
    //
    // #[test]
    // /// Normal run
    // fn run3(){
    //     init();
    //     info!("Test Nr. 2 -- Removing shit");
    //     let mut intervals: Vec<(u32, u32)> = vec![(1,20), (10,30), (8,15)];
    //     sort_vector(&mut intervals);
    //     let mut network = create_network_hashmap(& intervals);
    //     make_nested(&intervals, & mut network);
    //
    //     info!("run3 {:?}", network);
    //     assert_eq!(intervals.len(), 3);
    // }
    //
    // #[test]
    // /// Normal run
    // fn run4(){
    //     init();
    //     info!("Test Nr. 2 -- Removing shit");
    //     let mut intervals: Vec<(u32, u32)> = vec![(1,20), (10,30), (8,25), (0,40)];
    //     sort_vector(&mut intervals);
    //     let mut network = create_network_hashmap(& intervals);
    //     make_nested(&intervals, & mut network);
    //
    //     info!("run4 {:?}", intervals);
    //     assert_eq!(intervals.len(), 4);
    // }
    //

}
