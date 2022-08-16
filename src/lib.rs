use std::iter::FromIterator;
use hashbrown::{HashSet, HashMap};
use log::{debug, info, trace};
use smallvec::SmallVec;


/// Gives parent and child relation ship
#[derive(Eq, PartialEq, Debug, Clone)]
pub struct Network{
    pub parent: HashSet<(u32, u32)>,
    pub child: HashSet<(u32, u32)>,
}


impl Network{
    /// Create empty network
    pub fn new() -> Self{
        let parent: HashSet<(u32, u32)> = HashSet::new();
        let children: HashSet<(u32, u32)> = HashSet::new();
        Self{
            parent: parent,
            child: children,
        }
    }

}


/// Get all parents of a index
pub fn get_parents(start: &(u32, u32), hm: &HashMap<(u32, u32), Network>) -> HashSet<(u32, u32)>{
    let mut hh: HashSet<(u32, u32)> = HashSet::new();
    if hm.get(start).unwrap().parent.len() != 0{
        let j: HashSet<(u32, u32)> = hm.get(start).unwrap().parent.iter().cloned().collect();
        hh.extend(& j);
        for x in hm.get(start).unwrap().parent.iter(){
            let o = get_parents(x, hm);
            hh.extend(& o);
        }
    }
    hh
}

// /// Faster way to make the nestedness
// ///
// /// Still need parets
// pub fn make_nested_simple(intervals_sorted: & Vec<(u32, u32)>, order: & mut HashMap<(u32, u32), Network>){
//     let mut open_interval: (u32, u32) = (0,0);
//     for (start, end) in intervals_sorted.iter(){
//         if open_interval == (0,0){
//             open_interval = (start.clone(), end.clone());
//
//         }
//         else {
//             let result = fast_helper(&open_interval, &(start.clone(), end.clone()), order);
//             if result.len() == 1{
//                 order.get_mut(&(start.clone(), end.clone())).unwrap().parent.insert(result[0]);
//                 order.get_mut(&result[0]).unwrap().child.insert((start.clone(), end.clone()));
//             }
//             open_interval = (start.clone(), end.clone());
//         }
//     }
// }


// /// Recuive function for
// pub fn fast_helper(old: &(u32, u32), new: &(u32, u32), hm: & mut HashMap<(u32, u32), Network>) -> Vec<(u32, u32)>{
//     let mut ss: Vec<(u32, u32)> = Vec::new();
//     if (new.0 >= old.0) & (new.1 <= old.1){
//         ss.push(old.clone());
//     } else if hm.get(old).unwrap().parent.len() == 0 {
//         info!("nothing")
//     } else {
//         let parent = hm.get(old).unwrap().parent[0].clone();
//         ss.append(& mut fast_helper(& parent, new, hm))
//
//     }
//     ss
// }


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
pub fn make_nested(intervals_sorted: & Vec<(u32, u32)>, order: & mut HashMap<(u32, u32), Network>){
    trace!("Running nested");
    let mut open_intervals = SmallVec::<[(u32, u32); 20]>::new();
    
    // Iterate over sorted unique interval vector 
    for (start, end) in intervals_sorted.iter(){


        // If empty (just at the start), add to openintervals
        if open_intervals.len() == 0{
            open_intervals.push((start.clone(), end.clone()));
        } else {
            // Hits = direct relationship (child)
            // Overlaps = Overlap but not fully covered
            let mut hits: HashSet<(u32, u32)> = HashSet::new();
            let mut overlaps: HashSet<(u32, u32)> = HashSet::new();

            // Iterate over open intervals
            // For each interval: Start checking, when and if hit or overlap
            for (oldstart, oldend) in open_intervals.iter(){
                let mut hs = HashSet::new();
                let solution = checker_rec2(&(oldstart.clone(), oldend.clone()), &(start.clone(), end.clone()), order, false, & mut hs);
                hits.extend(& solution.0);
                overlaps.extend(& solution.1);

                //info!("HITS\t{:?}", hits);
                //info!("OVERLAPS\t{:?}", overlaps);
            }
            if hits.len() == 0{
                //info!("Nothing happening");
            }
            else {
                //info!("We need to filter");

                trace!("Overlaps {:?}", overlaps);
                trace!("Hits {:?}", hits);
                trace!("Start + End {:?}", (start, end));
                filter_hit(& mut hits);
                for x in hits{
                    order.get_mut(&(start.clone(), end.clone())).unwrap().parent.insert(x);
                    order.get_mut(&x).unwrap().child.insert((start.clone(), end.clone()));
                }


            }

            open_intervals =  SmallVec::<[(u32, u32); 20]>::new();
            for x in overlaps.iter() {
                open_intervals.push(x.clone());
            }
            open_intervals.push((start.clone(), end.clone()));

        }
       //info!("OP {:?}", &open_intervals);
        //info!("{:?}", order);

    }

    //info!("{:?}", order);
}


/// In a list of intervals, check if some are related to each other (very slow)
/// Complexity: O(n^2)
/// Input:
///     - candidates Vec<(u32, u32)>: List of intervals to check
///
/// Mut:
///     - candidates Vec<(u32, u32)>: List of intervals to check
///
/// Output:
///     --> Chaneging candicates
///
/// Running:
///     - Iterate over the list, check if the relation between the candidates
///
/// TODO:
/// - Stuff is removed at the end, checking like everything
/// PROBLEM:
/// - O(n^2) complexity
pub fn filter_hit(candidates: &mut HashSet<(u32, u32)>) {
    trace!("Running filter hit - Number of candidates {}", candidates.len());
    let cand_vector: Vec<(u32, u32)> = candidates.iter().cloned().collect();

    let mut remove_hashset = HashSet::new();

    for (i1, x) in cand_vector.iter().enumerate(){
        for y in cand_vector[i1+1..].iter(){
            // x is außerhalb von y (oder andersrum)
            // THIS IS NOT TRUE LOL
            if (x.0 <= y.0) & (x.1 >= y.1) {
                remove_hashset.insert(x.clone());
                // x is der parent
            } else if (x.0 >= y.0) & (x.1 <= y.1){
                remove_hashset.insert(y.clone());
            }

        }
    }
    //info!("Remove {:?}", remove_list);
    for x in remove_hashset.iter(){
        //trace!("tt {}", x-i);
        candidates.remove(x);
    }
    trace!("Running filter hit - Number of candidates {}", candidates.len());

}


/// In a list of intervals, check if some are related to each other (very slow)
/// Complexity: O(n^2)
/// Input:
///     - candidates Vec<(u32, u32)>: List of intervals to check
///
/// Mut:
///     - candidates Vec<(u32, u32)>: List of intervals to check
///
/// Output:
///     -
///
/// Running:
///     - Iterate over the list, check if the relation between the candidates
///
/// TODO:
/// - Stuff is removed at the end, checking like everything
/// - if true --> remove rerun
pub fn filter_hit2(candicates: &mut Vec<(u32, u32)>) {
    debug!("Running filter hit");
    let mut cand_index: Vec<usize> = (0..candicates.len()).collect();
    let mut detete = SmallVec::<[usize; 10]>::new();
    loop {
        let mut trigger = false;
        let mut remove_list = SmallVec::<[usize; 10]>::new();
        for i1 in cand_index.iter() {
            let n1 = candicates[*i1];
            for i2 in cand_index[i1 + 1..].iter() {
                let n2 = candicates[*i2];
                if (n1.0 <= n2.0) == (n1.1 >= n2.1) {
                    trigger = true;
                    if (n1.0 <= n2.0) & (n1.1 >= n1.0) {
                        remove_list.push(i1.clone());
                    } else {
                        remove_list.push(i2.clone());
                    }
                }
                if trigger {
                    break;
                }
            }

            if trigger {
                break;
            }
            remove_list.push(i1.clone());

        }

        for x in remove_list.iter(){
            cand_index.remove(x.clone().clone());
            detete.push(x.clone().clone());
        }
        if !trigger{
            break
        }
    }
    detete.sort();
    for (i, x) in detete.iter().enumerate(){
        candicates.remove(x - i);
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
pub fn checker_rec2(old: &(u32, u32), new: &(u32, u32), hm: & mut HashMap<(u32, u32), Network>, overlaps_parent: bool, hs: & mut HashSet<(u32, u32)>) -> (HashSet<(u32, u32)>, HashSet<(u32, u32)>) {
    debug!("Running checker recusive");
    //info!("Checking this interval {:?}", old);
    let mut hits: HashSet<(u32, u32)> = HashSet::new();
    let mut overlaps: HashSet<(u32, u32)> = HashSet::new();

    let mut now_overlapping: bool = false;

    // It is a hit!
    if (old.0 <= new.0) & (old.1 >= new.1) {

        hits.insert((old.0.clone(), old.1.clone()));
    } else {
        if (!overlaps_parent) & (old.1 > new.0) {
            overlaps.insert((old.0.clone(), old.1.clone()));
            now_overlapping = true;
        }
        if hm.get(old).unwrap().parent.len() != 0 {
            let mut vecc_p = Vec::new();
            for x in hm.get(old).unwrap().parent.iter() {
                //u.append(&mut helper1(x, new, hm));
                vecc_p.push(x.clone());
            }
            for x in vecc_p.iter() {
                if !hs.contains(x){
                    hs.insert(x.clone());
                    //info!("{:?}", x);
                    let jo = &mut checker_rec2(x, new, hm, now_overlapping, hs);
                    //info!("{:?}", jo);
                    hits.extend(& jo.0.clone());
                    overlaps.extend(& jo.1.clone());
                }
            }
        }
    }
    (hits.clone(), overlaps.clone())
}






/// Interval vector sorting
pub fn sort_vector(intervals: &mut Vec<(u32, u32)>){
    intervals.sort_by(|a, b| (a.0.cmp(&b.0).then(b.1.cmp(&a.1))));

}


/// Create the "network" for parent-child relationship
pub fn create_network_hashmap(intervals_sorted: & Vec<(u32, u32)>) -> HashMap<(u32, u32), Network>{
    let mut order: HashMap<(u32, u32), Network> = HashMap::new();
    for (start, end) in intervals_sorted.iter(){
        order.insert((start.clone(), end.clone()), Network::new());
    }
    order
}



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
    use hashbrown::HashSet;
    use log::info;
    use crate::{sort_vector, make_nested, create_network_hashmap, remove_duplicates, filter_hit};


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
        info!("{:?}", network)


    }
    #[test]
    fn filter1(){
        init();
        info!("Test Nr. 2 -- Removing shit");
        let intervals: Vec<(u32, u32)> = vec![(1,10), (2,4), (7,9), (8,9), (1,20)];
        let mut ii: HashSet<(u32, u32)> = intervals.iter().cloned().collect();
        filter_hit(& mut ii);

        info!("{:?}", ii);
        assert_eq!(ii.len(), 2);
    }

    #[test]
    fn filter2(){
        init();
        info!("Test Nr. 2 -- Removing shit");
        let intervals: Vec<(u32, u32)> = vec![(1,10), (2,4), (7,9), (1,2), (1,20)];
        let mut ii: HashSet<(u32, u32)> = intervals.iter().cloned().collect();
        filter_hit(& mut ii);

        info!("{:?}", ii);
        assert_eq!(ii.len(), 3);
    }


}
