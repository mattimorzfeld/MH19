function RunEDATuning(Nes,SetUps)

inflAll = .9:.1:2;
for jj=1:length(Nes)%[20 50 100 200]
    Ne = Nes(jj);
    for kk=1:length(SetUps)
        SetUp = SetUps(kk);
        for Gap = [2 4 6 8 10 12 14]
            TuneEDA(Ne,inflAll,SetUp,Gap);
        end
    end
end    