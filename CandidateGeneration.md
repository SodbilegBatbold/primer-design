# Introduction #

We start by declaring an object of the **primer\_pair** class, which we call **myprimers**:

```
  primer_pair myprimers;
```
Then we need to set the location and length of the primers. There are two ways to this. The first way is to set the location of the target sequence and the upstream/downstream regions in which we will search for primers. Also set the length range for the primers using the following:
```
  myprimers.set_target_location(200, 800);
  myprimers.set_flank_lengths(200, 200);
  myprimers.set_primer_length_range(18, 23);
```

Alternatively, we can set the location and lengths for the forward and reverse primers independently.
```
  myprimers.forward_primer.set_primer_location_range(0, 200);
  myprimers.reverse_primer.set_primer_location_range(800, 1000); 
  myprimers.forward_primer.set_primer_length_range(18, 23);
  myprimers.reverse_primer.set_primer_length_range(18, 23);
```

We can then generate a predefined number (presently 100) of candidate primers subject to hard constraints using: -
```
  myprimers.generate_candidates(template_sequence);
```

This is the minimum required to generate primer candidates and we can now analyse and sort these candidates to find the most suitable primer pair.