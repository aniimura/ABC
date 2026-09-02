import fs from 'fs';
const p = 'D:/Math_ABC3/ResearchPaper/blocked-leaves.json';
const j = JSON.parse(fs.readFileSync(p, 'utf8'));
const key = '[GenEll] Lemma 3.5(Vélu の商の半安定性——判別式の恒等式と p ∣ l の形式群)';
const e = j.blocked.find(b => b.key === key);
if (!e) { console.log('not found'); process.exit(1); }
e.hcopDoesNotDescend2026_09_02 =
  '★★☆**重要な訂正（2026-09-02、第 1408）——`hcop` は底変換で保たれない**。'
  + '☆`notDvdLClosed2026_09_02` で「残る仮定 `hcop`（`l ∤ jExp p E`）は消費側の '
  + '`PrimeToLocalHeights l` から来る」と書いたが、'
  + '★`VeluQuotOK` は `E.fld` の**任意の有限拡大 `M`** の素点 `P` で半安定性を要求し、'
  + '底変換では `jExp P (E⊗M) = e(P∣p)·jExp p E` となるので、'
  + '`l ∤ jExp p E` から `l ∤ jExp P (E⊗M)` は**出ない**（`l ∣ e(P∣p)` なら壊れる）。'
  + '☆したがって悪い素点でも **`l ∣ v_P(q)` の場合を正面から扱う**必要がある。'
  + '★その場合の核は `μ_l` 型でなく `⟨ζ^a q^{1/l}⟩` 型であり、'
  + '商は `E_{q′}`（`q′ = ζ^a q^{1/l}`、`v(q′) = v(q)/l ∈ ℤ`）でやはり乗法還元になる。'
  + '☆要るのは `veluQuotientFull (tateCurveAt q) ⟨q^{1/l}⟩` の Tate の計算'
  + '（第 1057 `veluQuotientFull_tate_mu` の「もう一方の部分群」版）である。';
e.remaining2026_09_02_final =
  '★★★★**2026-09-02 末の残り（第 1408 時点）**。'
  + '☆`p ∤ l` かつ `l ∤ jExp p E` の場合は**完全に閉じた**'
  + '（`semistableAt_veluQuot_of_not_dvd`、`Found/GenEll/VeluNotDvdL.lean`）。'
  + '★残るのは 3 つ:'
  + '（a）**`l ∣ v_P(q)` の悪い素点**——核が `μ_l` 型でない。'
  + '`veluQuotientFull (tateCurveAt q) ⟨q^{1/l}⟩` の Tate の計算が要る（5-10 ブロック）。'
  + '（b）**`p ∣ l`**——形式群か、`pDivLDualIdea2026_09_02` の双対同種の道（5-15 ブロック）。'
  + '（c）**インタフェース**——上の 2 つが閉じれば `VeluQuotOK` は仮定なしで出る（2-4 ブロック）。';
fs.writeFileSync(p, JSON.stringify(j, null, 2) + '\n', 'utf8');
console.log('recorded blocked-leaves');
