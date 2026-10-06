%% ========================================================================
%  bl.m - Analisi di stabilita' di un DAG (Tangle IOTA)
%  Autore: Andrea Martilotti Matricola 797428
%
%  Pipeline:
%   1. Parsing CSV e Hash Mapping -> matrice di adiacenza sparsa A
%   2. PageRank (metodo delle potenze, correzione dangling + teletrasporto)
%   3. Centralita' di Katz su una griglia di alpha
%   4. Concordanza dei nodi dominanti Katz-PageRank (top-k)
%   5. Laplaciana normalizzata simmetrica e autovettore di Fiedler
%   6. Sweep Cut incrementale e minimizzazione della conduttanza
%   7. Plateau del vettore di Fiedler e sua conduttanza
%   8. Grafici, salvataggio risultati, riepilogo CSV
%
%  Tutto cio' che lo script afferma e' garantito da un risultato teorico:
%   - Katz: soluzione esatta del sistema (nilpotenza -> serie di Neumann finita)
%   - PageRank: convergenza del metodo delle potenze (Perron-Frobenius, errore ~ alpha^k)
%   - lambda2 = 0  <=>  grafo non connesso
%   - Cheeger: lambda2/2 <= Phi_G <= Phi_sweep <= sqrt(2*lambda2), con lo sweep
%     eseguito su D^-1/2 v2 e su tutti i tagli (kMin = 1)
%  Non viene effettuata alcuna classificazione dell'attacco.
%
%  Convenzioni:
%   A(i,j) = 1 se la transazione i approva la transazione j (archi verso il passato)
%   verso = 1 : G = A'  (genesi -> tip, direzione dei random walker MCMC)
%   verso = 2 : G = A   (tip -> genesi)
%   Katz:  x = beta*(I - alpha*G')^(-1)*1
%          con verso = 2 coincide con l'eq. 5.1 (accumulo del cumulative weight)
%
%  Requisiti: MATLAB R2021a o successivo, nessun toolbox aggiuntivo.
% =========================================================================
clear; close all; clc;
rng(0);                                   % layout dei grafici riproducibile

%% 1. PARAMETRI -----------------------------------------------------------
n        = 500;                           % prime n transazioni del CSV
verso    = 1;                             % 1 = in avanti, 2 = all'indietro
dataset  = "dataset/healty/dataset_dag.csv";
%"C:\Users\andre\OneDrive\Desktop\dag_finale\dataset\healty\dataset_dag.csv"
cartella = "pdf";
colonnaTempo = "";                        % nome della colonna timestamp nel CSV ("" = ordine delle righe).
                                          % Se indicata, le righe vengono ordinate per tempo prima di
                                          % prendere le prime n: l'indice del nodo diventa il rango temporale.
maxNodiPlot=20000
% Katz
alphaGrid  = 1:-0.001:0;
beta       = 1;
nDominanti = 4;
% PageRank
alphaPR    = 0.85;
tolPR      = 1e-8;
maxIterPR  = 10000;
% Sweep Cut
kMin       = 1;
nBinPlateau = 200;

if verso == 1, pref = "av"; else, pref = "ind"; end
if ~isfolder(cartella), mkdir(cartella); end
diary off;
[percorsoDs, nomeDs] = fileparts(dataset);
[~, tipoDs] = fileparts(percorsoDs);
tag = tipoDs + "_" + nomeDs;
fileLog = fullfile(cartella, "log_exec_" + tag + ".txt");
if isfile(fileLog), delete(fileLog); end
diary(fileLog);
fprintf('Dataset: %s | n = %d | verso = %d\n', dataset, n, verso);
tempi = struct();

%% 2. CARICAMENTO E COSTRUZIONE DEL DAG ------------------------------------
tic;
colonneHash = {'Hash_Tx', 'Approvazione_1_Hash', 'Approvazione_2_Hash'};
opts = detectImportOptions(dataset);
opts = setvartype(opts, colonneHash, 'string');
if strlength(colonnaTempo) > 0
    opts.DataLines = [2, Inf];
    t = readtable(dataset, opts);
    assert(ismember(colonnaTempo, string(t.Properties.VariableNames)), ...
        'Colonna temporale "%s" non trovata nel CSV.', colonnaTempo);
    t = sortrows(t, colonnaTempo);
    t = t(1:min(n, height(t)), :);
else
    opts.DataLines = [2, n + 1];
    t = readtable(dataset, opts);
end

hashUnici = unique(t.Hash_Tx, 'stable');  % Hash Mapping: hash -> indice del nodo
N = numel(hashUnici);
if N < height(t)
    warning('%d righe con hash duplicato ignorate.', height(t) - N);
end
[~, locA]  = ismember(t.Hash_Tx,             hashUnici);
[~, locC1] = ismember(t.Approvazione_1_Hash, hashUnici);
[~, locC2] = ismember(t.Approvazione_2_Hash, hashUnici);

righe   = [locA;  locA];
colonne = [locC1; locC2];
validi  = colonne > 0 & righe ~= colonne;
A = spones(sparse(righe(validi), colonne(validi), 1, N, N));% binaria: doppie approvazioni contate una volta

if verso == 1, G = A'; else, G = A; end
nArchi = nnz(G);
Gd = digraph(G);
assert(isdag(Gd), 'Il grafo contiene cicli: non e'' un DAG.');
ordTopo = toposort(Gd);
prof = profonditaDAG(Gd, ordTopo);
GdAlt = digraph(A');%archi approvata -> approvante
altezza = altezzeDAG(GdAlt, toposort(GdAlt));%h(v): cammino piu' lungo dalla genesi
genesi  = find(full(sum(A, 2)) == 0);%transazioni che non approvano nulla
fprintf('Nodi: %d | Archi: %d | Profondita'': %d | Indice di nilpotenza: %d\n', ...
    N, nArchi, prof, prof + 1);
tempi.caricamento = toc;

%% 3. PAGERANK -------------------------------------------------------------
tic;
uno      = ones(N, 1);
outdeg   = full(sum(G, 2));
dangling = (outdeg == 0);
HT = (spdiags(1 ./ max(outdeg, 1), 0, N, N) * G)';   % trasposta della hyperlink matrix

% Prodotto con la Matrice di Google senza costruirla (resta O(|E|)):
% G_google' x = alpha*H' x + (alpha * massa_dangling + 1 - alpha)/N * 1
pr = uno / N;
convPR = false;
for iterPR = 1:maxIterPR
    xn = alphaPR * (HT * pr) + (alphaPR * sum(pr(dangling)) + 1 - alphaPR) / N;
    xn = xn / sum(xn);%corregge l'errore di arrotondamento
    delta = norm(xn - pr, 1);
    pr = xn;
    if delta < tolPR, convPR = true; break; end
end
if ~convPR
    warning('PageRank non converge in %d iterazioni (delta = %.2e).', maxIterPR, delta);
end
[valDom, dominanti] = maxk(pr, nDominanti);
fprintf('PageRank: convergenza in %d iterazioni (delta = %.2e)\n', iterPR, delta);
fprintf('Nodi dominanti PageRank: %s\n', mat2str(dominanti(:).'));
fprintf('Punteggi: %s\n', mat2str(valDom(:).', 4));
tempi.pagerank = toc;

%% 4. CENTRALITA' DI KATZ AL VARIARE DI ALPHA -----------------------------
% Iterazione x_{k+1} = alpha*G'*x_k + beta*1 con x_0 = beta*1 (serie di Neumann).
% G' e' nilpotente di indice prof+1: dopo prof passi la soluzione e' ESATTA
% per ogni alpha, senza criterio d'arresto. Il riscalamento congiunto di x e
% beta evita l'overflow senza cambiare la soluzione.
% Pari merito: per alpha sotto la soglia di crescita dei cammini (~1/2) i punteggi
% saturano a beta/(1-2*alpha) e le differenze ~(2*alpha)^h scendono sotto eps:
% migliaia di nodi risultano uguali in virgola mobile e maxk sceglie gli indici piu'
% bassi. In quella zona il ranking non e' determinato e viene escluso.
tic;
GT = G';
tolPari    = 1e-12;

Nalpha     = numel(alphaGrid);
indexKatz  = zeros(Nalpha, nDominanti);
nPari      = zeros(Nalpha, 1);
residuoMax = 0;
for kk = 1:Nalpha
    a = alphaGrid(kk);
    x = katzNeumann(a, beta, GT, prof, N);
    x = x / sum(x);
    residuoMax = max(residuoMax, residuoKatz(x, a, GT));
    [~, idx] = maxk(x, nDominanti);
    indexKatz(kk, :) = idx(:).';
    nPari(kk) = nnz(x >= max(x) * (1 - tolPari));
end
katzOK = nPari <= nDominanti;
if any(katzOK)
    alphaSoglia = min(alphaGrid(katzOK));
else
    alphaSoglia = NaN;
end
fprintf('Katz: residuo relativo massimo sulla griglia = %.2e\n', residuoMax);
fprintf('Katz: ranking univoco per %d valori di alpha su %d (alpha >= %.3f)\n', ...
    nnz(katzOK), Nalpha, alphaSoglia);
tempi.katz = toc;

%% 5. CONCORDANZA KATZ-PAGERANK 
% Calcolata solo sugli alpha con ranking di Katz univoco
top = dominanti(:).';
if any(katzOK)
    concOrdinata = mean(all(indexKatz(katzOK, :) == top, 2)) * 100;
    concInsieme  = mean(all(sort(indexKatz(katzOK, :), 2) == sort(top, 2), 2)) * 100;
else
    concOrdinata = NaN;  concInsieme = NaN;
end
fprintf('Concordanza top-%d (stesso ordine):  %.2f %%\n', nDominanti, concOrdinata);
fprintf('Concordanza top-%d (stesso insieme): %.2f %%\n', nDominanti, concInsieme);

%% 6. LAPLACIANA NORMALIZZATA SIMMETRICA E VETTORE DI FIEDLER -------------
tic;
Gsym  = spones(G + G');% proiezione non orientata (rete P2P)
comp  = conncomp(graph(Gsym));
nComp = max(comp);
[~, big] = max(accumarray(comp(:), 1));
idxLCC   = find(comp(:) == big);
fuoriLCC = find(comp(:) ~= big);
fprintf('Componenti connesse: %d | nodi fuori dalla componente principale: %d\n', nComp, numel(fuoriLCC));
if nComp > 1
    fprintf('Grafo NON connesso: lambda2(G) = 0, Phi_G = 0. Fiedler calcolato sulla componente principale.\n');
end
Gs = Gsym(idxLCC, idxLCC);
nS = numel(idxLCC);
assert(nS >= 3, 'Componente connessa troppo piccola per l''analisi spettrale.');

dS     = full(sum(Gs, 2));
volTot = sum(dS);
Dm12   = spdiags(1 ./ sqrt(dS), 0, nS, nS);
Nadj   = Dm12 * Gs * Dm12;% D^-1/2 A D^-1/2
Nadj   = (Nadj + Nadj') / 2;

if nS <= 3000
    L = eye(nS) - full(Nadj);% L_sym = I - D^-1/2 A D^-1/2
    [V, Dv] = eig(L);
    [lam, o] = sort(diag(Dv));
    lambda2 = lam(2);  v2 = V(:, o(2));
    lambda3 = lam(3);
else
    % Shift-invert su L_sym + 1e-3*I: restituisce gli autovalori piu'
    % vicini a -1e-3, cioe' 0, lambda2, lambda3. 
    Lsym = speye(nS) - Nadj;
    [V, Dv] = eigs(Lsym, 3, -1e-3);
    [lam, o] = sort(diag(Dv));
    lambda2 = lam(2);  v2 = V(:, o(2));
    lambda3 = lam(3);
end
lambda2 = max(lambda2, 0);
if abs(lambda3 - lambda2) < 1e-8
    warning('lambda2 ha molteplicita'' > 1: il vettore di Fiedler non e'' univoco.');
end
y = Dm12 * v2;% D^-1/2 v2: vettore su cui si esegue lo sweep
tempi.fiedler = toc;

%% 7. SWEEP CUT INCREMENTALE ----------------------------------------------
% Ordinamento per VALORE di D^-1/2 v2, poi per ogni prefisso S_k:
%   |cut(S_k)| = vol(S_k) - 2*E(S_k)  -> tutti i tagli in O(|E|) totale
tic;
[ySorted, ord] = sort(y);
P      = Gs(ord, ord);
Ein    = cumsum(full(sum(triu(P, 1), 1)).');% archi interni a S_k
volS   = cumsum(dS(ord));
taglio = volS - 2 * Ein;
volMin = min(volS, volTot - volS);

kk = (1:nS).';
valido = (kk >= kMin) & (kk <= nS - kMin) & (volMin > 0);
phi = inf(nS, 1);
phi(valido) = taglio(valido) ./ volMin(valido);
[phiBest, kBest] = min(phi);

S1 = ord(1:kBest);  S2 = ord(kBest + 1:end);
if sum(dS(S1)) <= sum(dS(S2))% Gruppo1 = gruppo a volume minore (indipendente dal segno di v2)
    G1l = S1;  G2l = S2;
else
    G1l = S2;  G2l = S1;
end
Gruppo1 = idxLCC(G1l);
Gruppo2 = idxLCC(G2l);
n1 = numel(Gruppo1);  n2 = numel(Gruppo2);

% Verifica della disuguaglianza di Cheeger
tolC = 1e-9;
cheegerOK = (phiBest >= lambda2 / 2 - tolC) && (phiBest <= sqrt(2 * lambda2) + tolC);
if ~cheegerOK
    if kMin > 1
        warning('Cheeger non verificata: con kMin > 1 il limite superiore non e'' garantito.');
    else
        warning('Cheeger non verificata: controllare i dati di ingresso.');
    end
end
tempi.sweep = toc;

%% 8. RISULTATI DEL TAGLIO -------------------------------------------------
subTaglio    = Gsym(Gruppo1, Gruppo2);
nArchiTaglio = nnz(subTaglio);
assert(nArchiTaglio == taglio(kBest), 'Calcolo incrementale del taglio incoerente.');
[rr, cc] = find(subTaglio);
nBound1 = numel(unique(rr));
nBound2 = numel(unique(cc));

fprintf('RISULTATI SPETTRALI\n');
fprintf('lambda2 (connettivita'' algebrica) = %.6g\n', lambda2);
fprintf('Cheeger: %.4g <= Phi_G <= Phi_sweep = %.4g <= %.4g  -> %s\n', ...
    lambda2 / 2, phiBest, sqrt(2 * lambda2), string(cheegerOK));
fprintf('Taglio: k = %d su %d | Gruppo1 (volume minore): %d nodi | Gruppo2: %d nodi\n', kBest, nS, n1, n2);
fprintf('Archi di taglio: %d | nodi di confine: %d (G1), %d (G2)\n', nArchiTaglio, nBound1, nBound2);

% lambda2*L^2 e' confrontabile tra finestre di dimensione diversa
lambda2L2 = lambda2 * prof^2;
fprintf('lambda2 * L^2 = %.4g\n', lambda2L2);

%% 8b. PLATEAU DEL VETTORE DI FIEDLER --------------------------------------
% Gruppo di nodi con D^-1/2 v2 quasi costante (bin piu' popolato dell'istogramma):
% il modo spettrale lento li tratta come un blocco unico. Lo sweep non puo' isolarlo
% se sta all'interno dell'ordinamento, quindi la sua conduttanza viene calcolata a parte.
[cntY, bordiY] = histcounts(y, nBinPlateau);
[~, bP]  = max(cntY);
inPlat   = y >= bordiY(bP) & y <= bordiY(bP + 1);
plateau  = idxLCC(inPlat);
nPlat    = numel(plateau);
resto    = setdiff((1:N).', plateau);
archiPlat = nnz(Gsym(plateau, resto));
volPlat   = full(sum(sum(Gsym(plateau, :))));
phiPlat   = archiPlat / min(volPlat, full(sum(Gsym(:))) - volPlat);
domInPlat = ismember(dominanti(:), plateau);
fprintf('Plateau: %d nodi (atteso con distribuzione uniforme: %.0f) | indici %d-%d | altezza %d-%d\n', ...
    nPlat, nS / nBinPlateau, min(plateau), max(plateau), min(altezza(plateau)), max(altezza(plateau)));
fprintf('Plateau: archi verso il resto = %d | Phi(plateau) = %.4g\n', archiPlat, phiPlat);
fprintf('Dominanti PageRank nel plateau: %d su %d\n', nnz(domInPlat), nDominanti);

%% 9. GRAFICI -------------------------------------------------------------
blu   = [0.13 0.38 0.63];
rosso = [0.86 0.08 0.24];
verde = [0 0.7 0];

% 9.1 Nodi dominanti di Katz al variare di alpha
figure('Color', 'w', 'Position', [100, 100, 800, 500]); hold on;
colori = lines(nDominanti);
hK = gobjects(nDominanti, 1);
for i = 1:nDominanti
    hK(i) = plot(alphaGrid, indexKatz(:, i), '.', 'MarkerSize', 8, 'Color', colori(i, :));
end
if ~isnan(alphaSoglia)
    xline(alphaSoglia, '--k', sprintf('\\alpha = %.3f', alphaSoglia), 'LineWidth', 1.2);
end
xlabel('Fattore di attenuazione \alpha', 'FontWeight', 'bold');
ylabel('Indice del nodo', 'FontWeight', 'bold');
title('Nodi dominanti (Katz) al variare di \alpha', 'FontSize', 14);
legend(hK, "Posizione " + string(1:nDominanti), 'Location', 'best');
stileAssi(); xlim([-0.01, 1.01]); hold off;
esporta(fullfile(cartella, pref + "_katz_" + tag + ".pdf"));

% % 9.2 Distribuzione del PageRank per altezza topologica
% figure('Color', 'w', 'Position', [100, 100, 800, 450]); hold on;
% plot(altezza, pr, '.', 'LineWidth', 0.8, 'Color', blu);
% plot(altezza(dominanti), valDom, 'o', 'MarkerSize', 7, 'MarkerFaceColor', rosso, 'MarkerEdgeColor', rosso);
% xlabel('Altezza topologica del nodo', 'FontWeight', 'bold');
% ylabel('Punteggio PageRank', 'FontWeight', 'bold');
% title('Distribuzione del punteggio PageRank', 'FontSize', 14);
% legend({'Punteggio PageRank', 'Nodi dominanti'}, 'Location', 'best');
% stileAssi(); xlim([-0.5, max(altezza) + 0.5]); hold off;
% esporta(fullfile(cartella, pref + "_page_rank_" + dataset + ".pdf"));
% 10.3 Distribuzione del PageRank sull'indice temporale
figure('Color', 'w', 'Position', [100, 100, 800, 450]); hold on;
plot(1:N, pr, '-', 'LineWidth', 0.9, 'Color', blu);
plot(dominanti, valDom, 'o', 'MarkerSize', 7, 'MarkerFaceColor', rosso, 'MarkerEdgeColor', rosso);
xlabel('Indice temporale del nodo', 'FontWeight', 'bold');
ylabel('Punteggio PageRank', 'FontWeight', 'bold');
title('Distribuzione del punteggio PageRank', 'FontSize', 14);
legend({'Punteggio PageRank', 'Nodi dominanti'}, 'Location', 'best');
stileAssi(); xlim([0, N + 1]); hold off;
esporta(fullfile(cartella, pref + "_page_rank_" + tag + ".pdf"));
% 9.3 DAG con i due gruppi dello Sweep Cut su fasce separate (asse X: ordine temporale)
if N <= maxNodiPlot
    figure('Color', 'w', 'Position', [100, 100, 900, 500]);
    p = plot(Gd, 'Marker', 's', 'MarkerSize', 5, 'LineWidth', 1.2, 'EdgeColor', [0.8 0.8 0.8], ...
        'NodeColor', [0.7 0.7 0.7]);      % grigio: nodi fuori dalla componente principale
    title("Sweep Cut: gruppo a volume minore (rosso) e resto della rete (blu) - " + tag, 'Interpreter', 'none');
    yc = zeros(1, N);
    yc(Gruppo2) = rand(1, n2);            % fascia bassa
    yc(Gruppo1) = rand(1, n1) + 2;        % fascia alta
    p.XData = 1:N;                        % ordine temporale
    p.YData = yc;
    highlight(p, Gruppo2, 'NodeColor', blu,   'Marker', 'o');
    highlight(p, Gruppo1, 'NodeColor', rosso, 'Marker', 's', 'MarkerSize', 7);
    [s, d] = find(G);                     % archi che attraversano il taglio
    attraversa = (ismember(s, Gruppo1) & ismember(d, Gruppo2)) | ...
                 (ismember(s, Gruppo2) & ismember(d, Gruppo1));
    if any(attraversa)
        highlight(p, s(attraversa), d(attraversa), 'EdgeColor', verde, 'LineWidth', 2.5);
    end
    highlight(p, genesi, 'NodeColor', 'g', 'MarkerSize', 10, 'Marker', 'diamond');
    if N > 5000                           % troppi oggetti per un PDF vettoriale: raster ad alta risoluzione
        exportgraphics(gcf, fullfile(cartella, pref + "_DAG_" + tag + ".pdf"), ...
            'ContentType', 'image', 'Resolution', 300, 'BackgroundColor','black');
    else
        esporta(fullfile(cartella, pref + "_DAG_" + tag + ".pdf"));
    end
end

% 9.4 Vettore di sweep: distribuzione e profilo ordinato
figure('Color', 'w', 'Position', [150, 150, 1000, 450]);
subplot(1, 2, 1); hold on;
histogram(y(G2l), 30, 'FaceColor', blu,   'EdgeColor', 'w', 'FaceAlpha', 0.8);
histogram(y(G1l), 30, 'FaceColor', rosso, 'EdgeColor', 'w', 'FaceAlpha', 0.8);
title('Distribuzione di D^{-1/2} v_2', 'FontSize', 12);
xlabel('Valore della componente');  ylabel('Numero di nodi');
legend({'Gruppo2', 'Gruppo1 (volume minore)'}, 'Location', 'best');
stileAssi(); hold off;

subplot(1, 2, 2); hold on;
inG1 = ismember(ord, G1l);
plot(find(~inG1), ySorted(~inG1), '.', 'MarkerSize', 12, 'Color', blu);
plot(find(inG1),  ySorted(inG1),  '.', 'MarkerSize', 12, 'Color', rosso);
inPlatOrd = inPlat(ord);
plot(find(inPlatOrd), ySorted(inPlatOrd), '.', 'MarkerSize', 6, 'Color', verde);
xline(kBest + 0.5, '--k', 'Taglio Sweep Cut', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5);
title('Profilo ordinato (Sweep Cut)', 'FontSize', 12);
xlabel('Posizione nell''ordinamento');  ylabel('Valore');
stileAssi(); hold off;
sgtitle(sprintf('\\lambda_2 = %.3g   |   \\Phi_{sweep} = %.3g', lambda2, phiBest), ...
    'FontWeight', 'bold', 'FontSize', 14);
esporta(fullfile(cartella, pref + "_fiedler_" + tag + ".pdf"));

%% 10. SALVATAGGIO RISULTATI ----------------------------------------------
fprintf('Tempi [s]: caricamento %.3f | PageRank %.3f | Katz %.3f | Fiedler %.3f | Sweep %.3f\n', ...
    tempi.caricamento, tempi.pagerank, tempi.katz, tempi.fiedler, tempi.sweep);

save(fullfile(cartella, pref + "_risultati_" + tag + ".mat"), ...
    'alphaGrid', 'indexKatz', 'nPari', 'katzOK', 'pr', 'dominanti', 'lambda2', ...
    'phiBest', 'altezza', 'Gruppo1', 'Gruppo2', 'y', 'plateau', 'tempi');

riepilogo = table(string(dataset), verso, N, nArchi, prof, nComp, numel(fuoriLCC), ...
    lambda2, lambda2L2, lambda2 / 2, phiBest, sqrt(2 * lambda2), cheegerOK, n1, n2, nArchiTaglio, ...
    nPlat, archiPlat, phiPlat, alphaSoglia, concOrdinata, concInsieme, ...
    'VariableNames', {'dataset', 'verso', 'nodi', 'archi', 'profondita', 'componenti', ...
    'fuori_componente', 'lambda2', 'lambda2_L2', 'cheeger_inf', 'phi_sweep', 'cheeger_sup', 'cheeger_ok', ...
    'n_gruppo1', 'n_gruppo2', 'archi_taglio', 'n_plateau', 'archi_plateau', 'phi_plateau', ...
    'alpha_soglia_katz', 'concordanza_ordine', 'concordanza_insieme'});
fileRiep = fullfile(cartella, "riepilogo.csv");
if isfile(fileRiep)
    writetable(riepilogo, fileRiep, 'WriteMode', 'append');
else
    writetable(riepilogo, fileRiep);
end
diary off;

%% FUNZIONI LOCALI --------------------------------------------------------
function d = profonditaDAG(Gd, ordTopo)
    % Lunghezza (in archi) del cammino piu' lungo: A^(d+1) = 0.
    d = max(altezzeDAG(Gd, ordTopo));
end

function h = altezzeDAG(Gd, ordTopo)
    % h(v) = lunghezza del cammino piu' lungo che termina in v (programmazione dinamica sull'ordinamento topologico, O(|V| + |E|)).
    h = zeros(numnodes(Gd), 1);
    for v = ordTopo(:).'
        succ = successors(Gd, v);
        if ~isempty(succ)
            h(succ) = max(h(succ), h(v) + 1);
        end
    end
end

function x = katzNeumann(a, beta, GT, prof, N)
    % Serie di Neumann finita: x = beta * sum_{k=0}^{prof} (a*G')^k * 1
    b = beta;
    x = b * ones(N, 1);
    for k = 1:prof
        x = a * (GT * x) + b;
        m = max(x);
        if m > 1e200% riscalamento congiunto di x e del termine noto
            x = x / m;  b = b / m;
        end
    end
end

function r = residuoKatz(x, a, GT)
    % (I - a*G') x deve essere un vettore costante. Il residuo e' rapportato alla
    % scala di x e non a z: per alpha > 1/2 z e' minuscolo rispetto a x  e il rapporto con z misurerebbe solo l'arrotondamento.
    z = x - a * (GT * x);
    r = (max(z) - min(z)) / max(abs(x));
end

function stileAssi()
    grid on;
    ax = gca;  ax.GridLineStyle = ':';  ax.GridAlpha = 0.5;
    ax.Color='none';
    lg = findobj(gcf, 'Type', 'legend');
    set(lg, 'Color', 'none');
end

function esporta(percorso)
    exportgraphics(gcf, percorso, 'ContentType', 'vector', 'BackgroundColor','black');
end