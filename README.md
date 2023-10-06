# Nested_intervals 

Find all relations between the provided intervals/ranges. 

## Definition of relation:
- Range A can be inside range B (nested)
  - A is the child of B 
  - B is the parent of A

## The data
Ranges/Intervals are represented as a array of size 2 (start and end position). These intervals can be parity or fully overlap other intervals. 



## Question
How to find all relations between the intervals?

## Solution
1) The vector of intervals need to be sorted by start position (ascending) and end position (descending). 
2) We initialize a second vector called: open_intervals for each interval which "open" (not closed) yet. 
3) We iterate over the sorted vector of intervals.
    - We sort the open_intervals by start position (descending) and end position (ascending)
    - Check the closest interval which covers the new_array -> we have found the direct relation
    - Only possibly, another range can be parent is when this interval has smaller end position than "hit" and smaller start position than "hit" -> we have found the indirect relation
    - Adjust smallest
    - Child and parents are reported (child.id, parent.id)

Plot will follow. 