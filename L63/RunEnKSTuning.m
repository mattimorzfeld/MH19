function RunEnKSTuning(Ne)

inflAll = .9:.05:2;
for SetUp = 1:5
    for Gap = [2 4 6 8 10 12 14]
        TuneEnKS(Ne,inflAll,SetUp,Gap);
    end
end
