import fs from 'fs';
const p = 'D:/Math_ABC3/ResearchPaper/blocked-leaves.json';
const j = JSON.parse(fs.readFileSync(p, 'utf8'));
const key = '[GenEll] Lemma 3.5(Vélu の商の半安定性——判別式の恒等式と p ∣ l の形式群)';
const e = j.blocked.find(b => b.key === key);
if (!e) { console.log('not found'); process.exit(1); }
e.goodPrimeClosed2026_09_02 =
  '★★★★★★★★★★★★★★★★**良い素点（`p ∤ l`）も閉じた（2026-09-02、第 1403）**。'
  + '☆`Found/GenEll/VeluGoodPrime.lean` の `semistableAt_veluQuot_goodPrime`——'
  + '仮定は `SemistableAt p E`・`0 ≤ jExp p E`・`l` 奇素数・`p ∤ l` だけ。'
  + '★極小モデルへの正規化（`minDeltaExp = max(0,−jExp) = 0` ＋ `minDeltaExp_eq`）を'
  + '足しただけで、第 1402 の恒等式に繋がった。'
  + '☆残るのは（1）**悪い素点の局所パッケージの配管**（第 1388 は sorry 0 だが'
  + '分岐版・分裂の二者択一が要る——α の葉の第 1373-1383 と同型の作業）と'
  + '（2）**`p ∣ l`**（形式群）である。';
e.pDivLRevised2026_09_02 =
  '★★☆**`p ∣ l` の見通しを改めた（2026-09-02、第 1403）**——'
  + '「原文の条件 `[L:ℚ]+1 < l` から `k = 0`」では**足りない**。'
  + '☆`VeluQuotOK` は `E.fld` の**任意の有限拡大 `M`** で半安定性を要求するので、'
  + '`M` が激しく分岐すれば核の点は形式群に入りうる（`e = v_p(l)` が大きくなる）。'
  + '★恒等式から `v(Δ′_Vélu) = 12(l−1)k` が出るが、`v(j(E′)) ≥ 0`（良還元）と'
  + '突き合わせると `c₄′ = c₄ + 240v` に `v(c₄′) ≥ 4(l−1)k` という大きな相殺が要る。'
  + '☆その相殺はまさに形式群の現象である'
  + '（`x_{jP} ≈ j^{−2}z^{−2}` なので `Σ_j j^{−4} ≡ 0 (mod l)` で主要項が消える）。'
  + '★★★したがって `p ∣ l` は**形式群（または Néron モデル・有限平坦群スキーム）が本体**であり、'
  + '「鋭い捩れ限界」だけでは閉じない。';
fs.writeFileSync(p, JSON.stringify(j, null, 2) + '\n', 'utf8');
console.log('recorded blocked-leaves');
