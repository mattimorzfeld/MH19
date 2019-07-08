function RunEnKFsqTuning(Ne)

inflAll = .9:.05:2;
for SetUp = [2 4 5]
    for Gap = [2 4 6 8 10 12 14]
        TuneEnKFsq(Ne,inflAll,SetUp,Gap);
    end
end
