function [est_network,clust_mtx] = clusters_community_detection(conn_mtx)

conn_mtx = conn_mtx-diag(diag(conn_mtx));
[Ci, Q] = modularity_louvain_und_sign(conn_mtx);
for irep = 1:100
  [Ci_next, Q_next] = modularity_louvain_und_sign(conn_mtx);
  if Q_next > Q
    Q = Q_next; Ci = Ci_next;
  end
end
clust_mtx = Ci(:);
est_network = zeros(length(clust_mtx));
for m = 1:length(clust_mtx)-1
    for n = m+1:length(clust_mtx)
        if clust_mtx(m) == clust_mtx(n)
            est_network(m,n) = 1; est_network(n,m) = 1;
        end
    end
end
