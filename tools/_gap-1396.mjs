import fs from 'fs';
const p = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const j = JSON.parse(fs.readFileSync(p, 'utf8'));
j.formalGroupTorsionSharp20260902 =
  '★★★★★★★★★★★★★★★★**`p ∣ l` の良い素点に残るのは「形式群の鋭い捩れ限界」ただ 1 つ**'
  + '（2026-09-02、第 1393-1396 の測定）。'
  + '☆必要な形: `E` が `p` で良還元・整、`P` が位数 `l`（奇素数）で `v_p(x_P) = −2k < 0`（深さ `k`）ならば'
  + ' **`(l−1)·k ≤ v_p(l)`**。★これがあれば、原文の条件 `[L:ℚ] + 1 < l` と'
  + ' `v_p(l) ≤ e ≤ [L:ℚ]` から `k = 0`、すなわち**核に深い点は無い**が出る。'
  + '☆在庫にあるのは**粗い版**だけ: `valAtLeast_pointCoords_of_le`（第 1077、証明済み）が'
  + ' `k ≤ v_p(l)` を与える（`l²·x` が整であることから）。★これでは `k = 0` は出ない。'
  + '☆鋭い版の古典的な証明は 2 通り: (a) 形式群 `[l](z) = lz + … + (unit)z^{l^h}` の Newton 多角形'
  + '（Silverman IV.6.1／VII.3.4、mathlib に形式群の対数は無い）、'
  + '(b) `preΨ_l = l·∏(x − r)` の Newton 多角形——核の `(l−1)/2` 個の根が付値 `−2k`、'
  + '残りの根は整なので、`x^{d−m}` の係数の整性から `v_p(l) ≥ 2k·(l−1)/2 = (l−1)k`。'
  + '★(b) は分解体と付値の延長と Vieta が要る（mathlib に Newton 多角形は無い、2026-09-02 確認）。'
  + '★★★これが閉じれば `p ∣ l` は「核の座標が整」に落ち、'
  + '第 1385（`semistableAt_of_disc_pow_eq`）＋第 1386（`veluKernelNorm_mem_primeSubring`）で'
  + '`p ∤ l` と**同じ形**になる（残るのは判別式の恒等式 1 本だけ）。';
j.veluDiscIdentityRoutes20260902 =
  '★★★★★★★★★★★★**判別式の恒等式 `Δ(E)^l = Δ(E′)·N⁴` の道を 3 つ測った**'
  + '（2026-09-02、第 1396 時点）。☆どれも「新しい理論を 1 つ建てる」規模である。'
  + '（i）**解析（`ℂ`）**——一意化から格子曲線に落とす所までは在庫'
  + '（第 1390-1392、`disc_pow_eq_veluQuot_mul_lattice` 1 本に落ちている）。'
  + '★残りは `Δ(Λ)` の η 積公式と `℘′(w) = −σ(2w)/σ(w)⁴`——'
  + 'mathlib に Weierstrass σ 函数も Δ の q 展開も無い。'
  + '（ii）**代数（終結式）**——`N = ±Res(ψ_C, 4x³+b₂x²+2b₄x+b₆)`、'
  + '`Δ′ = Δ^l / Res^4`（Dewaghe）。証明には `X(x) − X(x′) = ∏_{P∈C}(x − x(x′+P))/ψ(x)²` と'
  + '`Δ = 16∏_{i<j}(e_i − e_j)²` を使う——分解体と因子の理論が要る。'
  + '（iii）**剰余体（標数 `p`）の Vélu**／**Néron–Ogg–Shafarevich**／**モジュラー多項式 `Φ_l`**'
  + '——いずれも mathlib に無い（2026-09-02 確認）。'
  + '★★★**(i) が最短**——`ℂ` 側の足場（第 1330-1335、`Uniformization.lean` 5700 行）が既にある。';
fs.writeFileSync(p, JSON.stringify(j, null, 2) + '\n', 'utf8');
console.log('recorded mathlib-gap');
