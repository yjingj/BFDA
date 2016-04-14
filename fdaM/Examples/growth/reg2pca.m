
%  initialize iterations

niter    = 4;
itercell = cell(niter,4);
nharm    = 1;
WfdPar   = fdPar(basisw, 2, 1);

N = 10;
xfd = D2hgtffd(1:N);
xmeanfd  = mean(xfd);
xmeanmat = eval_fd(agefine, xmeanfd)*ones(1,N);

nfine = length(agefine);
tmat = agefine*ones(1,N);

%  carry out the iterative registration

for iter=1:niter
    
    %  principal components analysis
    
    pcastr = pca_fd(xfd, nharm);
    
    %  construct model matrix
    
    xharmmat = eval_fd(agefine, pcastr.harmfd);
    xscores  = pcastr.harmscr;
    xhatmat  = xharmmat*xscores' + xmeanmat;
    xhatfd   = smooth_basis(agefine, xhatmat, hgtfdPar);
    
    %  register to the model
    
    [yfd, warpfd] = registerfd(xhatfd, xfd, WfdPar);
    
    %  compute deformations
    
    deformmat = eval_fd(agefine, warpfd) - tmat;
    
    %  save results
    
    itercell{iter, 1} = yfd;
    itercell{iter, 2} = deformmat;
    itercell{iter, 3} = xharmmat;
    itercell{iter, 4} = pcastr.varprop;
    
    %  plot results
    
    figure(1)
    subplot(2,1,1)
    plot(xfd)
    xlabel('')
    ylabel('\fontsize{13} xfd')
    subplot(2,1,2)
    plot(yfd)
    xlabel('')
    ylabel('\fontsize{13} yfd')
    
    figure(2)
    plot(agefine, deformmat, '-', [1,18], [0,0], 'r:')
    
    figure(3) 
    plot(agefine, xharmmat, '-', [1,18], [0,0], 'r:');
    
    pause;
    
    %  set up for the next iteration
    
    xfd      = yfd;
    xmeanfd  = mean(xfd);
    xmeanmat = eval_fd(agefine, xmeanfd)*ones(1,N);
    
end


%  plot final results

%  xfd and yfd's

figure(1)
subplot(1,1,1)
xmat = eval_fd(agefine, D2hgtffd(1:N));
yarray = zeros(nfine,N,niter);
for iter=1:niter
    yarray(:,:,iter) = eval_fd(agefine, itercell{iter,1});
end
for i=1:N
    plot(agefine, xmat(:,i), 'g-', [1,18], [0,0], 'r:')
    hold on
    for iter=1:niter
        plot(agefine, yarray(:,i,iter), '-')
    end
    hold off
    axis([1,18,-4,2])
    xlabel('\fontsize{13} Age (years)')
    ylabel('\fontsize{13} Acceleration (cm/yr/yr)')
    title(['\fontsize{16} Case ',num2str(i)])
    pause
end

%  harmonics

figure(2)
for j=1:nharm
    plot(agefine, itercell{1,3}(:,j), 'g-', [1,18], [0,0], 'r:')
    hold on
    for iter=2:niter
        plot(agefine, itercell{iter,3}(:,j), '-')
    end
    hold off
    xlabel('\fontsize{13} Age (years)')
    ylabel('\fontsize{13} Harmonic')
    title(['\fontsize{16} Harmonic ',num2str(j)])
    pause
end

%  deformation functions

figure(3)
for i=1:N
    deformtotal = zeros(nfine,1);
    for iter=1:niter
        deformtotal = deformtotal + itercell{iter,2}(:,i);
    end
    plot(agefine, deformtotal, 'g-', [1,18], [0,0], 'r:')
    hold on
    for iter=1:niter
        plot(agefine, itercell{iter,2}(:,i), '-')
    end
    hold off
    xlabel('\fontsize{13} Age (years)')
    ylabel('\fontsize{13} Deformation h(t) - t (years)')
    title(['\fontsize{16} Case ',num2str(i)])
    pause
end

%  variance accounted for

figure(4)
varpropmat = zeros(niter,nharm);
for iter=1:niter
    varpropmat(iter,:) = itercell{iter,4}';
end
plot(1:niter, varpropmat, 'o-')
axis([1,niter,0,1])






