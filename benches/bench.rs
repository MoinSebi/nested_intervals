use criterion::{black_box, criterion_group, criterion_main, Criterion};
use related_intervals::{make_nested, sort_vector};

use std::cmp::Ordering;


/// Create a vector of tuples
/// These are in order
fn create_vector_linear() -> Vec<(u32, u32)>{
    let mut vector = Vec::new();
    //vector.push((0, 1000000));
    for value in 1..10000{
        vector.push((value +1, value +2))
    }
    return vector

}


/// Create a vector of tuples
/// These are in order
fn create_vector_jumps() -> Vec<[u32; 3]>{
    let mut vector = Vec::new();
    //vector.push((0, 1000000));
    for value in 1..100000{
        if value % 5000 == 0{

            vector.push([value +1, value +200000, 1]);
        } else if value % 200 == 0{
            vector.push([value +1, value +50, 1]);

        }

        else {
            vector.push([value +1, value +2, 1]);

        }
    }
    sort_vector(&mut vector);
    return vector

}







/// Testing the stuff
fn test_nested(intervals: &Vec<[u32; 3]>){

    //sort_vector(&mut intervals);
    let mut network = make_nested(&intervals);
}




fn criterion_benchmark(c: &mut Criterion) {
    let mut intervals = create_vector_jumps();




    //c.bench_function("faster network", |b| b.iter(|| test()));
    c.bench_function("faster network", |b| b.iter(|| test_nested(&intervals)));
    // c.bench_function("faster network", |b| b.iter(|| testt2()));


}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);