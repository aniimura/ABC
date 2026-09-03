import fs from 'fs';
const p = 'D:/Math_ABC3/ResearchPaper/blocked-leaves.json';
const j = JSON.parse(fs.readFileSync(p, 'utf8'));
const key = '[GenEll] Lemma 3.5(Vélu の商の半安定性——判別式の恒等式と p ∣ l の形式群)';
const e = j.blocked.find(b => b.key === key);
if (!e) { console.log('not found'); process.exit(1); }
e.notDvdLClosed2026_09_02 =
  '★★★★★★★★★★★★★★★★★★★★**`p ∤ l` の側が完全に閉じた（2026-09-02、第 1403-1407）**。'
  + '☆`semistableAt_veluQuot_of_not_dvd`（`Found/GenEll/VeluNotDvdL.lean`）——'
  + '良い素点（第 1403）と悪い素点（第 1406）を `jExp p E` の符号で場合分けする。'
  + '★悪い素点の配管は（i）第 1404 で不分岐の仮定を `^e` に弱め、'
  + '（ii）第 1405 で分裂性を落とし（非分裂は不分岐 2 次拡大へ、α の葉の第 1383 と同型）、'
  + '（iii）第 1406 で `L_p(ζ_l)` の局所パッケージを `p` ごとに作った。'
  + '☆残る仮定は 2 つだけ:`hcop`（`l ∤ jExp p E`＝`PrimeToLocalHeights`、'
  + '消費側の `hdag_of_stableLine` は持っている）と '
  + '`hc4L`（`c₄(E/C) ≠ 0`、数学的には自動だがモジュラー多項式か同種不変性が要る）。';
e.pDivLDualIdea2026_09_02 =
  '★★☆**`p ∣ l` に形式群を使わない道の候補（2026-09-02、第 1407 で気付いた）**。'
  + '☆恒等式を**双対同種にも当てる**:`E″ = E′/φ(E[l])` は `[l]` の分解なので `E″ ≅ E`'
  + '（Vélu の正規化では `u = l`）であり、`v(Δ(E″)) = ±12·v_p(l) = ±12e`。'
  + '★恒等式 `Δ(E′)^l = Δ(E″)·N′⁴` に `v(Δ(E′)) = 12S`（`S = Σ k_P ≥ 0`）と'
  + '`v(N′) ≤ 0`（第 1396 の内訳:整な点で `0`、深い点で `−3k`）を入れると'
  + '`12lS = v(Δ(E″)) + 4v(N′) ≤ v(Δ(E″)) = ±12e`、'
  + 'すなわち符号が正なら **`l·S ≤ e = v_p(l)`** が出る。'
  + '☆これは「鋭い捩れ限界」`(l−1)k ≤ e` より強く、`S = (l−1)k` なら `l(l−1)k ≤ e`。'
  + '★要るのは（a）双対同種の同定（`E′/φ(E[l]) ≅ E`、Vélu の合成が `[l]`）と'
  + '（b）符号の決定である。☆**形式群は要らない**——恒等式（第 1402）と'
  + '`v(N) ≤ 0`（第 1394-1396）だけで回る可能性がある。★次のセッションで測ること。';
fs.writeFileSync(p, JSON.stringify(j, null, 2) + '\n', 'utf8');
console.log('recorded blocked-leaves');
