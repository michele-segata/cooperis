function [AMP, PHA]= CFG_mode_per_cell(amp, pha, MAX_PHASE, d_pha, d_amp, average, random, seedn)

K = numel(pha);
[M, N] = size(pha{1});

AMP=zeros(M, N);
PHA=zeros(M, N);

amp_edges=0:d_amp:1; %for discretization
pha_edges=0:d_pha:MAX_PHASE; %for discretization

for m=1:M
    for n=1:N
        amp_cells_at_IJ=zeros(1,K);
        pha_cells_at_IJ=zeros(1,K);
        for I=1:K
            amp_cells_at_IJ(I)=amp{I}(m,n);
            pha_cells_at_IJ(I)=pha{I}(m,n);
        end

        if average
            phase = mean(pha_cells_at_IJ);
        else
            phase = maxOccurringElement(pha_cells_at_IJ, random, seedn);
        end
        pop = closestIndex(pha_edges, phase);
        PHA(m,n)=pha_edges(pop);

        if average
            amplitude = mean(amp_cells_at_IJ);
        else
            amplitude = maxOccurringElement(amp_cells_at_IJ, random, seedn);
        end
        pop = closestIndex(amp_edges, amplitude);
        AMP(m,n)=amp_edges(pop);

    end
end
