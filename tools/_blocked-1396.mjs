import fs from 'fs';
const p = 'D:/Math_ABC3/ResearchPaper/blocked-leaves.json';
const j = JSON.parse(fs.readFileSync(p, 'utf8'));
const key = '[GenEll] Lemma 3.5(Vélu の商の半安定性——判別式の恒等式と p ∣ l の形式群)';
let e = j.blocked.find(b => b.key === key);
if (!e) { e = { key }; j.blocked.push(e); }
e.measuredAt = '2026-09-02（第 1393-1396）';
e.reason =
  '★`Skeleton/GenEll/VeluSemistable.lean` の `semistableAt_veluQuotientFull` は'
  + '`EllModuliWitness`（したがって `theorem_3_8`・`cor_4_3`・`cor_4_4`）に残るただ 1 つの葉である。'
  + '☆場合は 3 つに割れている:'
  + '（a）**悪い素点（`p ∤ l`）**——★**閉じた**（第 1388 `semistableAt_veluQuotient_bad`、sorry 0）。'
  + '（b）**良い素点（`p ∤ l`）**——★**判別式の恒等式 1 本**に落ちている'
  + '（第 1385 `semistableAt_of_disc_pow_eq` ＋ 第 1387）。'
  + '（c）**`p ∣ l`（良い素点だけ。悪い素点では第 1042 が `p ∤ l` を強制する）**'
  + '——★恒等式に加えて**形式群の鋭い捩れ限界**が要る。';
e.exactMissingLemma =
  '★★★**欠けているものを 2 本に特定した（2026-09-02、第 1396）**。'
  + '（1）**判別式の恒等式** `Δ(E)^l = Δ(E/C)·(∏_{P∈C∖{O}}(2y_P+a₁x_P+a₃))^4`。'
  + '☆`l = 3,5,7` で数値検証済み（`tools/velu-disc-check.py`）、'
  + '降下は済んでいて（第 1390-1392）**格子曲線の上の 1 本**'
  + '（`disc_pow_eq_veluQuot_mul_lattice`）だけが `sorry` である。'
  + '★中身は `Δ(Λ)` の η 積公式と `℘′(w) = −σ(2w)/σ(w)⁴`——'
  + 'mathlib に Weierstrass σ 函数も Δ の q 展開も無い（`mathlib-gap.json` の'
  + '`veluDiscIdentityRoutes20260902` に 3 つの道を測ってある）。'
  + '（2）**形式群の鋭い捩れ限界** `(l−1)·k ≤ v_p(l)`（深さ `k` の位数 `l` の点）。'
  + '☆在庫は粗い版 `k ≤ v_p(l)`（第 1077）だけで、原文の条件 `[L:ℚ]+1 < l` から'
  + '`k = 0` を出すには鋭い版が要る（`mathlib-gap.json` の'
  + '`formalGroupTorsionSharp20260902`）。';
e.progress2026_09_02 =
  '★★★★★★★★★★★★**（c）の側で進んだ（第 1393-1396）**。'
  + '☆第 1393 `dvd_valAdd_Delta_sub_minDeltaExp`：どのモデルでも `12 ∣ v_p(Δ) − minDeltaExp`。'
  + '☆第 1393 `dvd_minDeltaExp_of_disc_pow_eq`：`v_p(Δ(E)) = 0` と恒等式と `3 ∣ v_p(N)` から'
  + '`12 ∣ minDeltaExp p E′`。'
  + '☆第 1394-1396：**`3 ∣ v_p(N)` を無条件で得た**'
  + '（`three_dvd_valAdd_veluKernelNorm`、`Found/GaloisRep/VeluKernelNormVal.lean`）。'
  + '★道具は 2 つ:(i)`Ψ₂Sq` と `Φ₂` の Bézout 恒等式（係数は `Δ` を出す、第 1394）で'
  + '「良い素点では両方が同時に深くならない」、'
  + '(ii)**深い点は倍化で深さが変わらない**（第 1395）＋ Fermat の小定理で `2^{l−1}P = P`。'
  + '☆これで「整な点では `v_p(t_P) = 0`、深い点では `−3k`」が場合分けなしに言える。'
  + '★★★残るのは上の 2 本だけである。';
fs.writeFileSync(p, JSON.stringify(j, null, 2) + '\n', 'utf8');
console.log('recorded blocked-leaves');
