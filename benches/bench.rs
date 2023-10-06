use criterion::{black_box, criterion_group, criterion_main, Criterion};
use related_intervals::{make_nested2, make_nested22, sort_vector};

use std::cmp::Ordering;


/// Create a vector of tuples
/// These are in order
fn creat_hashset() -> Vec<(u32, u32)>{
    let mut vector = Vec::new();
    //vector.push((0, 1000000));
    for value in 1..10000{
        vector.push((value +1, value +2))
    }
    return vector

}


/// Create a vector of tuples
/// These are in order
fn creat_hashset2() -> Vec<[u32; 2]>{
    let mut vector = Vec::new();
    //vector.push((0, 1000000));
    for value in 1..10000{
        if value % 500 == 0{
            vector.push([value +1, value +2000]);
        } else if value % 200 == 0{
            vector.push([value +1, value +500]);

        }

        else {
            vector.push([value +1, value +2]);

        }
    }
    sort_vector2(&mut vector);
    return vector

}

/// Interval vector sorting
pub fn sort_vector2(intervals: &mut Vec<[u32; 2]>){
    intervals.sort_by(|a, b| (a[0].cmp(&b[0]).then(b[1].cmp(&a[1]))));

}






/// Run the new function
fn testt(){
    let mut intervals = creat_hashset2();
    //sort_vector(&mut intervals);
    let mut network = make_nested2(&intervals);
}

//
//
// /// Run the new function
// fn testt2(){
//     let mut intervals = creat_hashset2();
//     //sort_vector(&mut intervals);
//     let mut network = make_nested22(&intervals);
// }




fn criterion_benchmark(c: &mut Criterion) {




    //c.bench_function("faster network", |b| b.iter(|| test()));
    c.bench_function("faster network", |b| b.iter(|| testt()));
    // c.bench_function("faster network", |b| b.iter(|| testt2()));


}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);