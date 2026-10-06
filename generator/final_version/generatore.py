"""
=======================================================================
GENERATORE TANGLE IOTA v2 -- VERSIONE CORRETTA
=======================================================================

CORREZIONI RISPETTO ALLA v1:

  [BUG 1] Collasso in catena lineare / doppia conferma sullo stesso tip.
          CAUSA: in una simulazione puramente sequenziale ogni transazione
          consuma 2 tip e ne produce 1 -> il tip pool si svuota subito e
          resta un unico tip, referenziato due volte. Il grafo degenerava
          in una catena (blockchain lineare), non un Tangle.
          FIX: introdotto il MODELLO DI LATENZA DI POPOV (whitepaper IOTA).
          Le transazioni arrivano con processo di Poisson di rate LAMBDA;
          una transazione emessa al tempo t vede il Tangle SOLO com'era al
          tempo (t - H_LATENCY), a causa del ritardo di propagazione di
          rete. Piu' transazioni emesse nella stessa finestra di latenza
          selezionano quindi dallo stesso insieme di tip, generando tip
          CONCORRENTI. Numero stazionario di tip: L ~ 2*LAMBDA*H_LATENCY.
          NOTA: questo NON e' la "finestra scorrevole" vietata nella
          specifica (pescare a caso tra le ultime N transazioni). La
          latenza e' il meccanismo fisico che genera i tip; la selezione
          resta rigorosamente un random walk MCMC.

  [BUG 2] Underflow numerico nella softmax.
          CAUSA: exp(-alpha*(H_x - H_y)) con differenze di peso cumulativo
          nell'ordine delle centinaia collassa a 0.0 per TUTTI i candidati
          -> somma nulla -> ZeroDivisionError ai bivi reali (nella v1 non
          emergeva solo perche' non esistevano bivi).
          FIX: stabilizzazione log-sum-exp (sottrazione del massimo
          esponente prima dell'esponenziale). Matematicamente identica,
          numericamente stabile.

  [BUG 3] Conferme duplicate non controllate.
          FIX: vincolo esplicito Conferma_1_Hash != Conferma_2_Hash,
          verificato a posteriori su tutte le righe. Unica eccezione
          MATEMATICAMENTE INEVITABILE: la primissima transazione dopo la
          Genesi, quando nel grafo esiste un solo nodo referenziabile.

MODELLO IMPLEMENTATO:
  - Arrivi di Poisson (inter-arrivo esponenziale, rate LAMBDA).
  - Visibilita' ritardata: la tx emessa a t vede solo cio' che e' stato
    emesso entro (t - H_LATENCY).
  - Tip = transazione visibile senza approvatori visibili.
  - Tip Selection: due random walk MCMC INDIPENDENTI dalla Genesi, con
    P(x->y) = exp(-alpha*(H_x - H_y)) / sum_z exp(-alpha*(H_x - H_z))
    calcolata sugli approvatori diretti visibili di x. alpha = 0.15.
  - Cumulative Weight H_x aggiornato ricorsivamente fino alla Genesi.
  - Out-degree esatto = 2 (0 per la Genesi), aciclicita' strutturale.

MODALITA' (parametro MODE):
  "healthy"         -> Tangle sano di riferimento
  "parasite_chain"  -> sotto-tangle parassita quasi isolato (pochi bridge)
  "eclipse_split"   -> due comunita' bilanciate con bottleneck stretto
  "star_spam"       -> hub artificiali bersagliati da flooding di spam

Schema CSV di output (INVARIATO rispetto alla v1):
  Hash_Tx, Mittente, Destinatario, Importo,
  Conferma_1_Hash, Conferma_2_Hash, Timestamp
Nessuna colonna di etichetta: l'infezione e' rilevabile solo per via
strutturale (analisi spettrale / Fiedler).
"""

import csv
import hashlib
import math
import os
import random
from collections import defaultdict

# ----------------------------------------------------------------------
# CONFIGURAZIONE
# ----------------------------------------------------------------------
N_TRANSACTIONS = 500      # totale transazioni, Genesi inclusa
ALPHA = 0.15              # selettivita' del random walk MCMC
LAMBDA_RATE = 1.0         # rate di arrivo delle transazioni (tx per unita' di tempo)
H_LATENCY = 6.0           # ritardo di propagazione -> tip attesi ~ 2*LAMBDA*H = 12
SEED = 42

MODE = "healthy"          # "healthy" | "parasite_chain" | "eclipse_split" | "star_spam"

# --- parametri delle infezioni ---
PARASITE_SIZE_RATIO = 0.22     # quota di transazioni nel sotto-tangle parassita
PARASITE_EXTRA_BRIDGES = 0     # bridge aggiuntivi oltre ai 2 di ancoraggio
PARASITE_START_FRACTION = 0.50 # istante (frazione della timeline) di nascita del parassita

ECLIPSE_SPLIT_RATIO = 0.5      # bilanciamento tra le due comunita'
ECLIPSE_EXTRA_BRIDGES = 1      # bridge aggiuntivi oltre ai 2 di ancoraggio

STAR_HUB_SIZE = 4              # numero di hub artificiali (deve essere >= 2)
STAR_SPAM_RATIO = 0.35         # quota di transazioni di spam
STAR_SPAM_PROB = 0.90          # probabilita' che una tx spam colpisca due hub

BASE_TIMESTAMP = 1735689600    # 2025-01-01 00:00:00 UTC
TIME_SCALE = 5                 # secondi per unita' di tempo della simulazione

# ----------------------------------------------------------------------
# STATO GLOBALE DEL TANGLE
# ----------------------------------------------------------------------
approves = {}                    # tx -> (conf_1, conf_2)
approvers = defaultdict(list)    # tx -> [tx che la confermano]
cum_weight = {}                  # tx -> H_x
issue_time = {}                  # tx -> istante di emissione (tempo simulato)
order = []                       # ordine di emissione
wallets = [f"wallet_{i:04d}" for i in range(120)]


def reset_state():
    """Azzera lo stato globale (necessario per generare piu' dataset nello stesso processo)."""
    approves.clear()
    approvers.clear()
    cum_weight.clear()
    issue_time.clear()
    order.clear()


def sha256_hash(payload: str) -> str:
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


# ----------------------------------------------------------------------
# CUMULATIVE WEIGHT (ricorsivo esatto, invariato)
# ----------------------------------------------------------------------
def get_past_set(tx_hash):
    """Past-set completo di tx_hash (tx_hash inclusa), via DFS esaustiva sugli out-edge."""
    visited = set()
    stack = [tx_hash]
    while stack:
        node = stack.pop()
        if node is None or node in visited:
            continue
        visited.add(node)
        a, b = approves[node]
        if a is not None:
            stack.append(a)
        if b is not None:
            stack.append(b)
    return visited


def propagate_cumulative_weight(new_tx, conf_1, conf_2):
    """La nuova tx entra nel future-set di TUTTI gli antenati di conf_1 e conf_2: +1 a ciascuno."""
    cum_weight[new_tx] = 1
    ancestors = set()
    for c in (conf_1, conf_2):
        if c is not None:
            ancestors |= get_past_set(c)
    for a in ancestors:
        cum_weight[a] += 1


# ----------------------------------------------------------------------
# VISIBILITA' E TIP POOL (modello di latenza)
# ----------------------------------------------------------------------
def is_visible(tx, cutoff):
    return issue_time[tx] <= cutoff


def get_visible_tips(now, allowed=None):
    """
    Tip visibili all'istante `now`: transazioni emesse entro (now - H_LATENCY)
    che non risultano confermate da alcuna transazione a sua volta visibile.
    `allowed`, se fornito, confina il calcolo a un sottoinsieme di nodi
    (usato per far crescere sotto-componenti isolate nelle infezioni).
    """
    cutoff = now - H_LATENCY
    tips = []
    for tx in order:
        if not is_visible(tx, cutoff):
            continue
        if allowed is not None and tx not in allowed:
            continue
        confermata = False
        for a in approvers[tx]:
            if is_visible(a, cutoff) and (allowed is None or a in allowed):
                confermata = True
                break
        if not confermata:
            tips.append(tx)
    return tips


# ----------------------------------------------------------------------
# RANDOM WALK MCMC (con stabilizzazione log-sum-exp)
# ----------------------------------------------------------------------
def mcmc_random_walk(start, now, allowed=None, exclude_tip=None):
    """
    Un esploratore parte da `start` e cammina lungo gli in-edge (approvers)
    VISIBILI, scegliendo ad ogni bivio con probabilita'

        P(x -> y) = exp(-alpha*(H_x - H_y)) / sum_z exp(-alpha*(H_x - H_z))

    calcolata con log-sum-exp per evitare underflow. Il walk termina su un
    tip visibile. `exclude_tip`, se fornito, impedisce al walk di fermarsi
    su quel tip: se e' l'unico candidato terminale rimasto il walk lo
    accetta comunque (situazione poi gestita dal chiamante).
    """
    cutoff = now - H_LATENCY
    current = start
    for _ in range(10000):  # guardia anti-loop (il grafo e' aciclico, non dovrebbe servire)
        candidates = [
            y for y in approvers[current]
            if is_visible(y, cutoff) and (allowed is None or y in allowed)
        ]
        if not candidates:
            return current

        # Se procedere porterebbe solo al tip escluso, fermati prima.
        if exclude_tip is not None and len(candidates) == 1 and candidates[0] == exclude_tip:
            sotto = [
                z for z in approvers[exclude_tip]
                if is_visible(z, cutoff) and (allowed is None or z in allowed)
            ]
            if not sotto:
                return current

        h_x = cum_weight[current]
        log_w = [-ALPHA * (h_x - cum_weight[y]) for y in candidates]
        m = max(log_w)
        w = [math.exp(lw - m) for lw in log_w]   # log-sum-exp: stabile numericamente
        tot = sum(w)
        probs = [x / tot for x in w]
        current = random.choices(candidates, weights=probs, k=1)[0]
    return current


def select_two_distinct_tips(start, now, allowed=None, max_retry=25):
    """
    Due esploratori MCMC indipendenti dalla Genesi -> due tip DISTINTI.
    Se i due walk convergono sullo stesso tip si ritenta; in ultima istanza
    si rilancia il secondo walk con quel tip escluso. Restituisce None se
    non esistono almeno due tip distinti (il chiamante avanzera' il tempo).
    """
    tips_disponibili = get_visible_tips(now, allowed)
    if not tips_disponibili:
        return None

    if len(tips_disponibili) == 1:
        # Regime di bootstrap: unico tip visibile. La prima referenza e' il
        # tip, la seconda viene presa via walk troncato tra i nodi visibili.
        # E' il fan-out che fa crescere il tip pool verso il regime stazionario.
        tip_1 = tips_disponibili[0]
        tip_2 = random_visible_node(start, now, allowed, exclude=tip_1)
        if tip_2 is None or tip_2 == tip_1:
            return None
        return tip_1, tip_2

    tip_1 = mcmc_random_walk(start, now, allowed)

    for _ in range(max_retry):
        tip_2 = mcmc_random_walk(start, now, allowed)
        if tip_2 != tip_1:
            return tip_1, tip_2

    tip_2 = mcmc_random_walk(start, now, allowed, exclude_tip=tip_1)
    if tip_2 != tip_1:
        return tip_1, tip_2

    # Fallback finale, ancora pesato (mai uniforme): sceglie tra i tip
    # rimanenti con la stessa legge di decadimento esponenziale.
    altri = [t for t in tips_disponibili if t != tip_1]
    h_ref = max(cum_weight[t] for t in altri)
    log_w = [-ALPHA * (h_ref - cum_weight[t]) for t in altri]
    m = max(log_w)
    w = [math.exp(lw - m) for lw in log_w]
    tot = sum(w)
    tip_2 = random.choices(altri, weights=[x / tot for x in w], k=1)[0]
    return tip_1, tip_2


# ----------------------------------------------------------------------
# CREAZIONE TRANSAZIONI
# ----------------------------------------------------------------------
def create_genesis():
    g = sha256_hash(f"GENESIS-TANGLE-{MODE}")
    approves[g] = (None, None)
    cum_weight[g] = 1
    issue_time[g] = 0.0
    order.append(g)
    return g


def emit_transaction(conf_1, conf_2, t, idx):
    mittente = random.choice(wallets)
    destinatario = random.choice([w for w in wallets if w != mittente])
    importo = round(random.uniform(0.5, 10000.0), 6)

    payload = f"{mittente}|{destinatario}|{importo}|{conf_1}|{conf_2}|{t}|{idx}"
    tx = sha256_hash(payload)
    while tx in approves:
        payload += "|x"
        tx = sha256_hash(payload)

    approves[tx] = (conf_1, conf_2)
    approvers[conf_1].append(tx)
    if conf_2 != conf_1:
        approvers[conf_2].append(tx)
    issue_time[tx] = t
    order.append(tx)
    propagate_cumulative_weight(tx, conf_1, conf_2)

    return {
        "Hash_Tx": tx,
        "Mittente": mittente,
        "Destinatario": destinatario,
        "Importo": importo,
        "Approvazione_1_Hash": conf_1,
        "Approvazione_2_Hash": conf_2,
        "Timestamp": BASE_TIMESTAMP + int(t * TIME_SCALE),
    }


def next_arrival(t):
    """Inter-arrivo esponenziale (processo di Poisson di rate LAMBDA_RATE)."""
    return t + random.expovariate(LAMBDA_RATE)


def random_visible_node(start, now, allowed=None, exclude=None, p_stop=0.35):
    """
    Walk MCMC TRONCATO: percorre il grafo come il walk normale ma puo'
    fermarsi a una profondita' casuale, restituendo quindi anche un nodo
    interno (non necessariamente un tip).

    Serve SOLO in fase di bootstrap, quando il grafo e' cosi' giovane da
    esporre un unico tip: in quel regime la seconda referenza viene presa
    tra i nodi visibili gia' confermati. Il protocollo IOTA lo consente
    (le referenze non devono obbligatoriamente essere tip) ed e' proprio
    questo fan-out iniziale che moltiplica il tip pool, dopodiche' subentra
    il regime stazionario con selezione MCMC piena su tip.
    """
    cutoff = now - H_LATENCY
    current = start
    for _ in range(10000):
        candidates = [
            y for y in approvers[current]
            if is_visible(y, cutoff) and (allowed is None or y in allowed)
        ]
        if not candidates:
            break
        if random.random() < p_stop and current != start and current != exclude:
            break
        h_x = cum_weight[current]
        log_w = [-ALPHA * (h_x - cum_weight[y]) for y in candidates]
        m = max(log_w)
        w = [math.exp(lw - m) for lw in log_w]
        tot = sum(w)
        current = random.choices(candidates, weights=[x / tot for x in w], k=1)[0]

    if current == exclude or current == start:
        visibili = [
            x for x in order
            if is_visible(x, cutoff) and x != exclude
            and (allowed is None or x in allowed)
        ]
        if not visibili:
            return None
        h_ref = max(cum_weight[x] for x in visibili)
        log_w = [-ALPHA * (h_ref - cum_weight[x]) for x in visibili]
        m = max(log_w)
        w = [math.exp(lw - m) for lw in log_w]
        tot = sum(w)
        current = random.choices(visibili, weights=[x / tot for x in w], k=1)[0]

    return current


# ----------------------------------------------------------------------
# BOOTSTRAP: seme minimo che porta il tip pool a >= 2 tip
# ----------------------------------------------------------------------
def bootstrap(genesis, t, idx, rows):
    """
    La primissima transazione dopo la Genesi e' costretta a referenziare la
    Genesi due volte: nel grafo esiste un solo nodo, e' un vincolo
    matematico, non una scorciatoia. Da li' in avanti tutte le conferme
    sono distinte: la seconda transazione referenzia (tx1, Genesi), la
    terza (tx2, tx1), dopodiche' la latenza inizia a produrre tip
    concorrenti e subentra la selezione MCMC piena.
    """
    t = next_arrival(t)
    r1 = emit_transaction(genesis, genesis, t, idx)   # unico caso degenere ammesso
    rows.append(r1)
    idx += 1

    t = next_arrival(t)
    r2 = emit_transaction(r1["Hash_Tx"], genesis, t, idx)
    rows.append(r2)
    idx += 1

    t = next_arrival(t)
    r3 = emit_transaction(r2["Hash_Tx"], r1["Hash_Tx"], t, idx)
    rows.append(r3)
    idx += 1

    return t, idx


def genesis_row(genesis):
    return {
        "Hash_Tx": genesis,
        "Mittente": "GENESIS",
        "Destinatario": "GENESIS",
        "Importo": 0.0,
        "Approvazione_1_Hash": "",
        "Approvazione_2_Hash": "",
        "Timestamp": BASE_TIMESTAMP,
    }


# ----------------------------------------------------------------------
# MODE 1: TANGLE SANO
# ----------------------------------------------------------------------
def build_healthy(n):
    rows = []
    g = create_genesis()
    rows.append(genesis_row(g))

    t, idx = 0.0, 1
    t, idx = bootstrap(g, t, idx, rows)

    tip_counts = []
    while len(rows) < n:
        t = next_arrival(t)
        sel = select_two_distinct_tips(g, t)

        # FIX TEMPORALE: Nessun salto di netto. Se non ci sono tip, il tempo scorre piano.
        while sel is None:
            t += 0.1
            sel = select_two_distinct_tips(g, t)

        c1, c2 = sel
        tip_counts.append(len(get_visible_tips(t)))
        rows.append(emit_transaction(c1, c2, t, idx))
        idx += 1

    print(f"[INFO healthy] tip pool medio: {sum(tip_counts)/len(tip_counts):.2f} "
          f"(atteso ~ {2*LAMBDA_RATE*H_LATENCY:.1f}) | max: {max(tip_counts)}")
    return rows, g


# ----------------------------------------------------------------------
# MODE 2: PARASITE CHAIN
# ----------------------------------------------------------------------
def build_parasite_chain(n):
    """
    Il Tangle principale cresce sano; a meta' timeline nasce un sotto-tangle
    parassita ancorato al grafo principale da 2 sole transazioni di
    ancoraggio (+ eventuali bridge extra). Il parassita evolve poi con un
    proprio tip pool INTERNO, restando strutturalmente quasi isolato.
    """
    rows = []
    g = create_genesis()
    rows.append(genesis_row(g))

    t, idx = 0.0, 1
    t, idx = bootstrap(g, t, idx, rows)

    parasite_size = int(n * PARASITE_SIZE_RATIO)
    start_at = int(n * PARASITE_START_FRACTION)

    main_nodes = {r["Hash_Tx"] for r in rows}
    parasite_nodes = set()
    bridges = 0
    while len(rows) < n:
            t = next_arrival(t)

            # Controlliamo se siamo nel periodo in cui l'attaccante è attivo
            fase_parassita_attiva = (len(rows) >= start_at) and (len(parasite_nodes) < parasite_size)

            # Diamo all'attaccante una certa probabilità di emettere la transazione in questo istante
            # Ad esempio, l'attaccante genera il 25% delle transazioni in questa fase
            turno_attaccante = fase_parassita_attiva and (random.random() < 0.25)

            if not turno_attaccante:
                # --- CRESCITA DEL MAIN TANGLE (Continua sempre a scorrere) ---
                sel = select_two_distinct_tips(g, t, allowed=main_nodes)
                if sel is None:
                    t += H_LATENCY
                    continue
                c1, c2 = sel
                r = emit_transaction(c1, c2, t, idx)
                main_nodes.add(r["Hash_Tx"])
                rows.append(r)
                idx += 1
                continue

            # --- ANCORAGGIO DEL PARASSITA (Solo le prime 2 tx) ---
            if len(parasite_nodes) < 2:
                sel = select_two_distinct_tips(g, t, allowed=main_nodes)
                if sel is None:
                    t += H_LATENCY
                    continue
                c1, c2 = sel
                r = emit_transaction(c1, c2, t, idx)
                parasite_nodes.add(r["Hash_Tx"])
                rows.append(r)
                idx += 1
                bridges += 1
                continue

            # --- CRESCITA INTERNA DEL PARASSITA ---
            p_root = min(parasite_nodes, key=lambda x: issue_time[x])
            sel = select_two_distinct_tips(p_root, t, allowed=parasite_nodes)
            if sel is None:
                t += H_LATENCY
                continue
            c1, c2 = sel

            # Bridge occasionali (come nel tuo codice originale)
            if bridges < 2 + PARASITE_EXTRA_BRIDGES and random.random() < 0.04:
                main_tips = get_visible_tips(t, allowed=main_nodes)
                if main_tips:
                    c2 = random.choice(main_tips)
                    bridges += 1

            r = emit_transaction(c1, c2, t, idx)
            parasite_nodes.add(r["Hash_Tx"])
            rows.append(r)
            idx += 1

    print(f"[INFO parasite_chain] principale: {len(main_nodes)} | "
          f"parassita: {len(parasite_nodes)} | bridge: {bridges}")
    return rows, g


# ----------------------------------------------------------------------
# MODE 3: ECLIPSE SPLIT (CORRETTO - COMPATIBILE CON LA VERIFICA)
# ----------------------------------------------------------------------
def build_eclipse_split(n):
    rows = []
    g = create_genesis()
    rows.append(genesis_row(g))

    t, idx = 0.0, 1
    t, idx = bootstrap(g, t, idx, rows)

    comm_a = {r["Hash_Tx"] for r in rows}
    comm_b = set()

    # ========================================================
    # BOOTSTRAP DELLA COMUNITA' B (Senza conferme duplicate)
    # ========================================================
    t = next_arrival(t)
    sel = None
    while sel is None: # Aspettiamo che ci siano tip disponibili in A
        t += 1.0
        sel = select_two_distinct_tips(g, t, allowed=comm_a)
    anchor1, anchor2 = sel

    # 1. Il vero e proprio Bridge (Root di B) agganciato ad A
    r1_b = emit_transaction(anchor1, anchor2, t, idx)
    comm_b.add(r1_b["Hash_Tx"])
    rows.append(r1_b)
    idx += 1

    # 2. Secondo nodo di B (approva il Root e l'anchor di A) -> Niente duplicati!
    t = next_arrival(t)
    r2_b = emit_transaction(r1_b["Hash_Tx"], anchor1, t, idx)
    comm_b.add(r2_b["Hash_Tx"])
    rows.append(r2_b)
    idx += 1

    # 3. Terzo nodo di B (approva il secondo nodo e il Root) -> Niente duplicati!
    t = next_arrival(t)
    r3_b = emit_transaction(r2_b["Hash_Tx"], r1_b["Hash_Tx"], t, idx)
    comm_b.add(r3_b["Hash_Tx"])
    rows.append(r3_b)
    idx += 1

    # Avanziamo il tempo artificialmente per rendere visibili questi nodi
    t += H_LATENCY

    b_root = r1_b["Hash_Tx"]

    # ========================================================
    # CRESCITA PARALLELA REALE
    # ========================================================
    while len(rows) < n - ECLIPSE_EXTRA_BRIDGES:
        t = next_arrival(t)

        if random.random() < ECLIPSE_SPLIT_RATIO:
            # Turno della Comunità A
            sel = select_two_distinct_tips(g, t, allowed=comm_a)
            if sel is not None:
                c1, c2 = sel
                r = emit_transaction(c1, c2, t, idx)
                comm_a.add(r["Hash_Tx"])
                rows.append(r)
                idx += 1
        else:
            # Turno della Comunità B
            sel = select_two_distinct_tips(b_root, t, allowed=comm_b)
            if sel is not None:
                c1, c2 = sel
                r = emit_transaction(c1, c2, t, idx)
                comm_b.add(r["Hash_Tx"])
                rows.append(r)
                idx += 1

    # ========================================================
    # BRIDGE FINALI (RICUCITURA DELL'ECLIPSE)
    # ========================================================
    while len(rows) < n:
        t = next_arrival(t)
        tips_a = get_visible_tips(t, allowed=comm_a)
        tips_b = get_visible_tips(t, allowed=comm_b)

        if not tips_a or not tips_b:
            t += H_LATENCY
            continue

        c1 = mcmc_random_walk(g, t, allowed=comm_a)
        c2 = mcmc_random_walk(b_root, t, allowed=comm_b)
        if c1 == c2:
            continue

        r = emit_transaction(c1, c2, t, idx)
        comm_a.add(r["Hash_Tx"])
        rows.append(r)
        idx += 1

    print(f"[INFO eclipse_split] A: {len(comm_a)} | B: {len(comm_b)} | bridge: {2 + ECLIPSE_EXTRA_BRIDGES}")
    return rows, g


# ----------------------------------------------------------------------
# MODE 4: STAR SPAM
# ----------------------------------------------------------------------
def build_star_spam(n):
    rows = []
    g = create_genesis()
    rows.append(genesis_row(g))

    t, idx = 0.0, 1
    t, idx = bootstrap(g, t, idx, rows)

    nodes = {r["Hash_Tx"] for r in rows}
    start_at = int(n * 0.40) # L'attacco spam inizia al 40% del DAG

    # FASE INIZIALE SANA
    while len(rows) < start_at:
        t = next_arrival(t)
        sel = select_two_distinct_tips(g, t, allowed=nodes)
        while sel is None:
            t += 0.1
            sel = select_two_distinct_tips(g, t, allowed=nodes)
        c1, c2 = sel
        r = emit_transaction(c1, c2, t, idx)
        nodes.add(r["Hash_Tx"])
        rows.append(r)
        idx += 1

    # CREAZIONE HUB
    hubs = []
    for _ in range(STAR_HUB_SIZE):
        t = next_arrival(t)
        sel = select_two_distinct_tips(g, t, allowed=nodes)
        while sel is None:
            t += 0.1
            sel = select_two_distinct_tips(g, t, allowed=nodes)
        c1, c2 = sel
        r = emit_transaction(c1, c2, t, idx)
        nodes.add(r["Hash_Tx"])
        hubs.append(r["Hash_Tx"])
        rows.append(r)
        idx += 1

    # CRESCITA PARALLELA (Traffico onesto miscelato con lo Spam)
    n_spam = 0
    while len(rows) < n:
        t = next_arrival(t)

        # Turno spammer: probabilità del 35% di iniettare spam
        if random.random() < 0.35:
            c1, c2 = random.sample(hubs, 2)
            n_spam += 1
        else:
            sel = select_two_distinct_tips(g, t, allowed=nodes)
            while sel is None:
                t += 0.1
                sel = select_two_distinct_tips(g, t, allowed=nodes)
            c1, c2 = sel

        r = emit_transaction(c1, c2, t, idx)
        nodes.add(r["Hash_Tx"])
        rows.append(r)
        idx += 1

    print(f"[INFO star_spam] sane: {len(nodes) - n_spam} | hub: {STAR_HUB_SIZE} | spam: {n_spam}")
    return rows, g


# ----------------------------------------------------------------------
# VERIFICA
# ----------------------------------------------------------------------
def verify(rows, genesis):
    hashes = {r["Hash_Tx"] for r in rows}
    assert len(hashes) == len(rows), "Hash duplicati!"

    pos = {h: i for i, h in enumerate(order)}
    duplicati = []
    for r in rows:
        h = r["Hash_Tx"]
        c1, c2 = r["Approvazione_1_Hash"], r["Approvazione_2_Hash"]
        if h == genesis:
            assert c1 == "" and c2 == "", "La Genesi non deve avere conferme!"
            continue
        assert c1 and c2, "Out-degree diverso da 2!"
        assert pos[c1] < pos[h] and pos[c2] < pos[h], "Violazione di aciclicita'!"
        assert issue_time[c1] < issue_time[h] and issue_time[c2] < issue_time[h], \
            "Violazione dell'ordinamento temporale!"
        if c1 == c2:
            duplicati.append(h)

    n_dup = len(duplicati)
    in_deg = [len(approvers[h]) for h in hashes]

    print(f"[VERIFICA] transazioni          : {len(rows)}")
    print(f"[VERIFICA] hash univoci         : {len(hashes)}")
    print(f"[VERIFICA] aciclicita'           : OK")
    print(f"[VERIFICA] out-degree = 2        : OK")
    print(f"[VERIFICA] conferme duplicate   : {n_dup} "
          f"({'solo bootstrap, inevitabile' if n_dup <= 1 else 'ANOMALIA!'})")
    print(f"[VERIFICA] in-degree max/medio  : {max(in_deg)} / {sum(in_deg)/len(in_deg):.2f}")
    print(f"[VERIFICA] H_x Genesi           : {cum_weight[genesis]}")
    assert n_dup <= 1, "Troppe conferme duplicate: il Tangle sta collassando in catena!"


# ----------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------
def generate(mode, out_path):
    global MODE
    MODE = mode
    reset_state()
    random.seed(SEED)

    builders = {
        "healthy": build_healthy,
        "parasite_chain": build_parasite_chain,
        "eclipse_split": build_eclipse_split,
        "star_spam": build_star_spam,
    }
    rows, g = builders[mode](N_TRANSACTIONS)
    verify(rows, g)

    fields = ["Hash_Tx", "Mittente", "Destinatario", "Importo",
              "Approvazione_1_Hash", "Approvazione_2_Hash", "Timestamp"]
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"-> scritto {out_path} ({len(rows)} righe)\n")


if __name__ == "__main__":
    # Compatibile sia con esecuzione da file .py sia con Jupyter/Colab,
    # dove __file__ non e' definito.
    try:
        base = os.path.dirname(os.path.abspath(__file__))
    except NameError:
        base = os.getcwd()

    for m in ["healthy", "parasite_chain", "eclipse_split", "star_spam"]:
        name = "dataset_dag_mcmc_v2.csv" if m == "healthy" else f"dataset_dag_mcmc_v2_{m}.csv"
        print("=" * 62)
        print(f"MODE: {m}")
        print("=" * 62)
        generate(m, os.path.join(base, name))