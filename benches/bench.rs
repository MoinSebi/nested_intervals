use hashbrown::{HashSet};
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use related_intervals::{create_network_hashmap, filter_hit, get_parents, make_nested, sort_vector};



fn makehs() -> Vec<(u32, u32)>{
    let mut o = Vec::new();
    for x in 1..100{
        o.push((x+1, x+10))
    }
    let mut ii: HashSet<(u32, u32)> = o.iter().cloned().collect();
    return o

}

fn test(){
    let mut intervals = makehs();
    //sort_vector(&mut intervals);
    let mut network = create_network_hashmap(& intervals);

    make_nested(&intervals, & mut network);
}

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("faster network", |b| b.iter(|| test()));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);