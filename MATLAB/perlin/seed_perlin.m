function P = seed_perlin()
    P = [randperm(256)' randperm(256)' randperm(256)']-1;
end