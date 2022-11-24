clear all
load('b_struct_brie_drowning')

b_struct.slr_rate = [10e-3*ones(1,12) 2e-3];

out = {barrier_model(b_struct)};
