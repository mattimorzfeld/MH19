function RunJPPFTuning(Nes,SetUps)

inflAll = .05:.05:.8;
for jj=1:length(Nes)%[20 50 100 200]
    Ne = Nes(jj);
    for kk=1:length(SetUps)
        SetUp = SetUps(kk);
        for Gap = [2 4 6 8 10 12 14]
            TuneJPPF(Ne,inflAll,SetUp,Gap);
        end
    end
end    