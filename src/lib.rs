use std::iter::FromIterator;
use std::collections::{HashSet, HashMap};
use std::hash::Hash;
use log::debug;


/// Each interval has a parent and child
#[derive(Hash, Eq, PartialEq, Debug, Clone)]
pub struct Network{
    pub parent: Vec<(u32, u32)>,
    pub child: Vec<(u32, u32)>,
}


impl Network{
    /// Create empty network
    pub fn new() -> Self{
        let parent: Vec<(u32, u32)> = Vec::new();
        let children: Vec<(u32, u32)> = Vec::new();
        Self{
            parent: parent,
            child: children,
        }
    }

}

/// Get all parents
///
/// Not sure if cloned is the right way
pub fn get_parents(start: &(u32, u32), hm: &HashMap<(u32, u32), Network>) -> HashSet<(u32, u32)>{
    let mut hh: HashSet<(u32, u32)> = HashSet::new();
    if hm.get(start).unwrap().parent.len() != 0{
        let mut j: HashSet<(u32, u32)> = hm.get(start).unwrap().parent.iter().cloned().collect();
        hh = hh.union(& mut j).cloned().collect();
        for x in hm.get(start).unwrap().parent.iter(){
            let mut o = get_parents(x, hm);
            hh = hh.union(& mut o).cloned().collect();
        }
    }
    hh
}

/// Faster way to make the nestedness
///
/// Still need parets
pub fn make_nested_simple(intervals_sorted: & Vec<(u32, u32)>, order: & mut HashMap<(u32, u32), Network>){
    let mut open_interval: (u32, u32) = (0,0);
    for (start, end) in intervals_sorted.iter(){
        if open_interval == (0,0){
            open_interval = (start.clone(), end.clone());

        }
        else {
            let result = fast_helper(&open_interval, &(start.clone(), end.clone()), order);
            if result.len() == 1{
                order.get_mut(&(start.clone(), end.clone())).unwrap().parent.push(result[0]);
                order.get_mut(&result[0]).unwrap().child.push((start.clone(), end.clone()));
            }
            open_interval = (start.clone(), end.clone());
        }
    }
}


/// Recuive function for
pub fn fast_helper(old: &(u32, u32), new: &(u32, u32), hm: & mut HashMap<(u32, u32), Network>) -> Vec<(u32, u32)>{
    let mut ss: Vec<(u32, u32)> = Vec::new();
    if (new.0 >= old.0) & (new.1 <= old.1){
        ss.push(old.clone());
    } else if hm.get(old).unwrap().parent.len() == 0 {
        println!("nothing")
    } else {
        let parent = hm.get(old).unwrap().parent[0].clone();
        ss.append(& mut fast_helper(& parent, new, hm))

    }
    ss
}


/// Give sense to Network 
/// 
/// Adding parent-child relation ship
/// This is also a wrapper for recursive function
/// 
pub fn make_nested(intervals_sorted: & Vec<(u32, u32)>, order: & mut HashMap<(u32, u32), Network>){
    debug!("Running nested");
    let mut open_intervals: Vec<(u32, u32)> = Vec::new();
    
    // Iterate over sorted unique interval vector 
    for (start, end) in intervals_sorted.iter(){


        // If empty (just at the start), add to openintervals
        if open_intervals.len() == 0{
            open_intervals.push((start.clone(), end.clone()));
        } else {
            let mut hits: Vec<(u32, u32)> = Vec::new();
            let mut overlaps: Vec<(u32, u32)> = Vec::new();
            
            // Iterate over open intervals
            // For each interval: Start checking, when and if hit or overlap
            for (oldstart, oldend) in open_intervals.iter(){
                let mut solution = checker_rec(&(oldstart.clone(), oldend.clone()), &(start.clone(), end.clone()), order, false);
                hits.append(& mut solution.0);
                overlaps.append(& mut solution.1);

                //println!("HITS\t{:?}", hits);
                //println!("OVERLAPS\t{:?}", overlaps);
            }
            if hits.len() == 0{
                //println!("Nothing happening");
            }
            else if hits.len() == 1{
                //println!("Only one");
                order.get_mut(&(start.clone(), end.clone())).unwrap().parent.push(hits[0]);
                order.get_mut(&hits[0]).unwrap().child.push((start.clone(), end.clone()));
            } else {
                //println!("We need to filter");
                filter_hit(& mut hits);
                for x in hits{
                    order.get_mut(&(start.clone(), end.clone())).unwrap().parent.push(x);
                    order.get_mut(&x).unwrap().child.push((start.clone(), end.clone()));
                }


            }

            open_intervals =  Vec::new();
            for x in overlaps.iter() {
                open_intervals.push(x.clone());
            }
            open_intervals.push((start.clone(), end.clone()));

        }
       //println!("OP {:?}", &open_intervals);
        //println!("{:?}", order);

    }

    //println!("{:?}", order);
}


/// Remove all parents here - maybe slow
pub fn filter_hit(candicates: &mut Vec<(u32, u32)>) {
    loop{
        let mut trigger = false;
        let mut remove_list: HashSet<usize> = HashSet::new();
        for (i1, x) in candicates.iter().enumerate(){
            for (i2, y) in candicates[i1+1..].iter().enumerate(){
                // x is außerhalb von y (oder andersrum)
                if (x.0<= y.0) == (x.1 >= y.1){
                    trigger = true;
                    // x is der parent
                    if (x.0<= y.0) & (x.1 >= y.0){
                        remove_list.insert(i1);
                    // y ist der parent
                    } else {
                        remove_list.insert(i2);
                    }
                }
            }
            if trigger{
                break;
            }
        }
        //println!("Remove {:?}", remove_list);
        for (i,x) in remove_list.iter().enumerate(){
            candicates.remove(x-i);
        }
        if !trigger{
            break;
        }
    }

}

/// Recursive
pub fn checker_rec(old: &(u32, u32), new: &(u32, u32), hm: & mut HashMap<(u32, u32), Network>, overlaps_parent: bool) -> (Vec<(u32, u32)>, Vec<(u32, u32)>) {
    debug!("Running checker recusive");
    //println!("Checking this interval {:?}", old);
    let mut hits: Vec<(u32, u32)> = Vec::new();
    let mut overlaps: Vec<(u32, u32)> = Vec::new();

    let mut now_overlapping: bool = false;

    // It is a hit!
    if (old.0 <= new.0) & (old.1 >= new.1) {

        hits.push((old.0.clone(), old.1.clone()));
    } else {
        if (!overlaps_parent) & (old.1 > new.0) {
            overlaps.push((old.0.clone(), old.1.clone()));
            now_overlapping = true;
        }

        if hm.get(old).unwrap().parent.len() != 0 {
            let mut vecc_p = Vec::new();
            for x in hm.get(old).unwrap().parent.iter() {
                //u.append(&mut helper1(x, new, hm));
                vecc_p.push(x.clone());
            }
            for x in vecc_p.iter() {
                //println!("{:?}", x);
                let jo = &mut checker_rec(x, new, hm, now_overlapping);
                //println!("{:?}", jo);
                hits.append(&mut jo.0.clone());
                overlaps.append(&mut jo.1.clone());
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
        eprintln!("Removed duplicates");
        intervals.retain(|e| uniques.insert(*e));
    }

}


/// Check if start and stop in order
///
/// Return: bool
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
    println!("LENNNN {}", change_list.len());
    println!("LENNNN {:?}", change_list);
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
    use crate::{sort_vector, make_nested, create_network_hashmap, remove_duplicates, filter_hit, check_overlapping, make_nested_simple, get_parents, start_stop_check, start_stop_order};


    // cargo test -- --nocapture
    #[test]
    fn basic() {
        // We test remove and and general function
        println!("\nRunning test 1");
        //assert_eq!(2 + 2, 4);
        let i1: (u32, u32) = (1, 10);
        let i2: (u32, u32) = (2, 4);
        let i3: (u32, u32) = (7,9);
        let i4: (u32, u32) = (7,9);
        let i5: (u32, u32) = (8,9);

        let mut intervals: Vec<(u32, u32)> = Vec::new();


        intervals.push(i1);
        intervals.push(i2);
        intervals.push(i3);
        intervals.push(i4);
        intervals.push(i5);

        remove_duplicates(&mut intervals);
        assert_eq!(intervals.len(), 4);
        sort_vector(&mut intervals);
        let mut network = create_network_hashmap(& intervals);

        make_nested(&intervals, & mut network);
        println!("{:?}", network)


    }
    #[test]
    fn hit_remover_testing(){
        println!("Test Nr. 2 -- Removing shit");
        let i1: (u32, u32) = (1, 12);
        let i2: (u32, u32) = (10, 12);
        let i3: (u32, u32) = (1, 6);
        let i4: (u32, u32) = (1, 3);
        let i5: (u32, u32) = (1, 2);
        let mut k: Vec<(u32, u32)> = Vec::new();
        k.push(i1);
        k.push(i2);
        k.push(i3);
        k.push(i4);
        k.push(i5);
        println!("{:?}", k);
        filter_hit(& mut k);
        println!("{:?}", k);
        assert_eq!(k.len(), 2);
    }

    #[test]
    fn overlap_check(){
        let i1: (u32, u32) = (1, 20);
        let i2: (u32, u32) = (1, 10);
        let i3: (u32, u32) = (7, 20);
        let i4: (u32, u32) = (11, 13);
        let i5: (u32, u32) = (8,9);
        let mut k: Vec<(u32, u32)> = Vec::new();
        k.push(i1);
        k.push(i2);
        k.push(i3);
        k.push(i4);
        k.push(i5);

        sort_vector(&mut k);
        let mut network = create_network_hashmap(& k);
        println!("{:?}", k);
        //make_nested(&k, & mut network);
        make_nested(&k, & mut network);
        println!("{:?}", network);
        let g = get_parents(&i2, &network);
        println!("parents {:?}", g);
    }

    #[test]
    fn overlap_check2(){
        let i1: (u32, u32) = (1, 20);
        let i2: (u32, u32) = (1, 10);
        let i6: (u32, u32) = (2,6);
        let i4: (u32, u32) = (3,4);
        let i5: (u32, u32) = (11,12);
        let mut k: Vec<(u32, u32)> = Vec::new();
        k.push(i1);
        k.push(i2);
        k.push(i4);
        k.push(i5);
        k.push(i6);

        sort_vector(&mut k);
        let mut network = create_network_hashmap(& k);
        println!("{:?}", k);
        //make_nested(&k, & mut network);
        make_nested(&k, & mut network);
        println!("{:?}", network);
        let g = get_parents(&i2, &network);
        println!("parents {:?}", g);
    }


    #[test]
    fn check_fast(){
        let i1: (u32, u32) = (1, 20);
        let i2: (u32, u32) = (1, 10);
        let i3: (u32, u32) = (9, 20);
        let i4: (u32, u32) = (11, 13);
        let i5: (u32, u32) = (11,12);
        let i6: (u32, u32) = (11,11);
        let mut k: Vec<(u32, u32)> = Vec::new();
        k.push(i1);
        k.push(i2);
        k.push(i3);
        k.push(i4);
        k.push(i5);


        let mut k2: Vec<(u32, u32)> = Vec::new();
        k2.push(i1);
        k2.push(i2);
        k2.push(i4);
        k2.push(i5);
        k2.push(i6);




        let k10 = check_overlapping(& mut k);
        let k11 = check_overlapping(& mut k2);
        assert_eq!(k10, true);
        assert_eq!(k11, false);


        let mut network = create_network_hashmap(& k2);
        println!("{:?}", k2);
        //make_nested(&k, & mut network);
        make_nested_simple(&k2, & mut network);
        println!("{:?}", network);
    }

    #[test]
    fn check_order() {
        let i1: (u32, u32) = (1, 20);
        let i2: (u32, u32) = (10, 1);
        let i3: (u32, u32) = (9, 20);
        let i4: (u32, u32) = (11, 13);
        let i5: (u32, u32) = (11, 12);
        let i6: (u32, u32) = (11, 11);
        let mut k: Vec<(u32, u32)> = Vec::new();
        k.push(i1);
        k.push(i2);
        k.push(i3);
        k.push(i4);
        k.push(i5);

        let mut k2: Vec<(u32, u32)> = Vec::new();
        k2.push(i6);
        k2.push(i2);
        k2.push(i3);
        k2.push(i4);
        k2.push(i5);

        assert_eq!(start_stop_check(&k), false);
        assert_eq!(start_stop_check(&k2), false);
        start_stop_order(& mut k);
        assert_eq!(k[1], (1,10));

    }
}
