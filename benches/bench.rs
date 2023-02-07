use criterion::{black_box, criterion_group, criterion_main, Criterion};
use related_intervals::{create_network_hashmap, filter_hit, make_nested, make_nested2, sort_vector};

use std::cmp::Ordering;



fn creat_hashset() -> Vec<(u32, u32)>{
    let mut o = Vec::new();
    o.push((0,1000000));
    for x in 1..10000{
        o.push((x+1, x+2))
    }
    return o

}




fn test(){
    let mut intervals = creat_hashset();
    //sort_vector(&mut intervals);
    let mut network = create_network_hashmap(& intervals);

    make_nested(&intervals, & mut network);
}

fn testt(){
    let mut intervals = creat_hashset();
    //sort_vector(&mut intervals);
    let mut network = make_nested2(&intervals);
}




fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("faster network", |b| b.iter(|| test()));
    c.bench_function("faster network", |b| b.iter(|| testt()));


}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);