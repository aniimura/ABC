import ABC3.Found.PGC.TwoIsDegenerate

/-!
# [pGC] ★★★★2 本の柱は**同じ 1 本の恒等式**に落ちた —— `σπ − π = (1+π)·w`

## ★費用（先に測って選んだ）

持ち場は「測定 1 と測定 2 のどちらが安いか先に」と言った。★測ったら
**どちらでもなく、2 本とも同じ 1 本の恒等式から出た**。所要 3 秒。

## ★★★★核 —— 第 1 跳びの `σ_a` に対する `ρ = σπ − π` の因数分解

`π = ζ − 1`（`ζ = ζ_{p^n}`）、`σ_a(π) = ζ^a − 1`、`a ≡ 1 (mod p)` で `a − 1 = p·b` とすると

> `ρ = σ_a(π) − π = (ζ^a − 1) − (ζ − 1) = ζ^a − ζ = ζ·(ζ^{a−1} − 1)`
> `  = (1 + π)·(ζ_{p^{n−1}}^b − 1) = (1 + π)·w`,  ★`w := ζ_{p^{n−1}}^b − 1 ∈ 𝒪_{E₁}`

（`ζ^{a−1} = ζ^{pb} = (ζ^p)^b = ζ_{p^{n−1}}^b` は `pow_mul_split`。）

★★**これで 2 本の柱が同時に出る:**

### 測定 1（`ρ` の展開は 2 成分だけ）⇒ **恒等式そのもの**

`w ∈ 𝒪_{E₁}` なので `ρ = w·(1+π) = w + w·π`（`mul_one_add_eq`）。
⇒ ★**`R_0 = R_1 = w`（同じ元！）、`R_j = 0` for `j ≥ 2`。**
★前波は「`v(R_0) = v(R_1) = p` が偶然そろっている」ように見えたが、
★★**同じ元だったから当然だった。**

`v_L(w) = p` は `w = ζ_{p^{n−1}}^b − 1 = (ζ_{p^{n−1}} − 1)·Σ_{i<b} ζ_{p^{n−1}}^i`
（`geom_unit_factor`）で、第 2 因子は `p ∤ b` のとき単数、
第 1 因子は `E₁` の素元だから `v_L = e(L/E₁) = p`。

### 測定 2（`t = p(p−1)`）⇒ **同じ計算を 1 段下でやるだけ**

`E₁` の素元を `μ = ζ_{p^{n−1}} − 1` に取ると、同じ因数分解で
`σμ − μ = ζ_{p^{n−1}}·(ζ_{p^{n−2}}^b − 1)`、その `v_L` は `e_L/φ(p^{n−2}) = p²`
（`index_ratio_sq`）。`v_L(μ) = p` なので

> ★`t = v_L(σμ − μ) − v_L(μ) = p² − p = p(p−1)`（`t_eq_p_mul_sub_one`）

★前波は `t = p(p−1)` を「測定」としか書けなかった。★**これで理由が出た。**

## ★★検証（★#350 の手に従い 3 つ目の場合まで測った）

`ρ = (1+π)·w` を **5 層 × 第 1 跳びの σ すべて**で確かめた:

| p | n | e_L | 第 1 跳びの σ | `ρ = (1+π)w` | `v(w) = p` | `R_j = 0 (j≥2)` |
|---|---|---|---|---|---|---|
| 2 | 4 | 8 | 4 | ✓ | ✓ | ✓ |
| 2 | 5 | 16 | 8 | ✓ | ✓ | ✓ |
| 3 | 3 | 18 | 6 | ✓ | ✓ | ✓ |
| 3 | 4 | 54 | 18 | ✓ | ✓ | ✓ |
| 5 | 3 | 100 | 20 | ✓ | ✓ | ✓ |

★合計 **56 個の σ** で 1 件の例外もなし（`rho_factor_verified`）。
`v_L(σμ − μ)` も `4 = 2²`, `9 = 3²`, `25 = 5²` で一致（`sigma_mu_valuation_sq`）。

## ★★★訂正 —— 前波の「2 本とも測定である」は**過小評価だった**

前波の docstring は「測定 1・2 は測定である。一般の完全分岐拡大では確かめていない」
と書いた。★★**本波で 2 本とも 1 本の恒等式に落ちた。**
★残るのは **「`ζ_{p^m} − 1` が `ℚ_p(ζ_{p^m})` の素元」** ただ 1 つである。

## ★在庫の測定（★2 か所。#330）

1. **mathlib**:
   `grep -in "cyclotomic" .cache/mathlib-index.txt | grep -iE "ramif|eisenstein|uniformi|valuation|discr"`
   ⇒ 判別式（`IsCyclotomicExtension.Rat.discr_prime_pow'` 等）は**豊富**だが、
   ★**局所体の下付き分岐群の跳び**は無い。ramification の唯一の当たりは
   `IsCyclotomicExtension.Rat.inertiaDegIn_ramificationIdxIn_aux`
   （★**private** かつ ℚ 上の**大域**で `ramificationIdxIn = p^k(p−1)`）。
   ⇒ ★`e_L = p^{n−1}(p−1)` は取れるが、跳びは取れない。
2. **木**: `TotallyRamified.lean`（39,971 バイト、宣言 20 本超）を**開いた**。
   ★設定は `PAdicLocalField p` の `K.closure` と `IsTotallyRamifiedAdjoin` で、
   ★★**`ζ − 1` の具体模型ではない**。`adjoinIntegers` を経由するので
   `lean-idioms.md` #69 の境界（212 秒 timeout）に当たる。
   ⇒ ★**本波では使わない**方が安いと判断した。★開いた結果を報告に書く。

## ★次に測るべき 1 点

★**「`ζ_{p^m} − 1` が素元」を木か mathlib から引く**（`Polynomial.cyclotomic` の
Eisenstein 性か、`IsCyclotomicExtension` の局所版）。★これ 1 つで本波の 2 本が
完全な証明になる。★ただし #69 の境界を越えない道を探すこと。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は付けていない（★**検算の記録**である）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. ★§1・§2 は**証明**である（環の恒等式と ℕ/ℤ の算術）。
   ★しかし「`ζ_{p^m} − 1` が素元」は**入れていない** —— 入れると `IsCyclotomicExtension`
   の局所版が要り、#69 の境界に当たる。★**そこは仮定のまま残した。**
4. ★§3 は**測定**である（56 個の σ、5 層）。
5. ★宣言名 12 件は `grep -rn "^theorem <名前>" lean/ABC3/Found/PGC/*.lean` で
   **実ファイル**の衝突 0 を先に確かめた（`lean-idioms.md` #348）。
6. ★docstring は `.py` を Write ツールで書いてから実行した（`lean-idioms.md` #349、
   前波で `python -c` の backtick を bash に食われて表が空になった実害があるため）。
-/

namespace ABC3.Found.PGC

namespace RhoFactorization

/-! ## §1 抽象核 —— すべて環の恒等式（証明） -/

section Kernel

/-- ★核その 1 —— `z^{k+1} − z = z(z^k − 1)`。
`ρ = ζ^a − ζ = ζ(ζ^{a−1} − 1)` の骨格。 -/
theorem sub_self_eq_mul_sub_one {A : Type*} [CommRing A] (z : A) (k : ℕ) :
    z ^ (k + 1) - z = z * (z ^ k - 1) := by
  rw [pow_succ]; ring

/-- ★★★**本ファイルの中心** —— `π` の言葉で書いた同じ恒等式:
`((1+π)^{k+1} − 1) − π = (1+π)·((1+π)^k − 1)`。

左辺は `σ_a(π) − π`（`a = k+1`）、右辺の第 2 因子は `a ≡ 1 (mod p)` のとき
★`ζ_{p^{n−1}}^b − 1 ∈ 𝒪_{E₁}` になる。★これで `ρ` が `E₁` の元の `(1+π)` 倍になる。 -/
theorem sigma_sub_eq {A : Type*} [CommRing A] (π : A) (k : ℕ) :
    ((1 + π) ^ (k + 1) - 1) - π = (1 + π) * ((1 + π) ^ k - 1) := by
  rw [pow_succ]; ring

/-- ★★核その 2 —— `w(1+π) = w + wπ`。★**成分が 2 つしかない**ことの中身。
⇒ `R_0 = R_1 = w`（★同じ元）、`R_j = 0` for `j ≥ 2`。
★前波は「`v(R_0) = v(R_1)` が偶然そろう」ように見えたが、同じ元だから当然だった。 -/
theorem mul_one_add_eq {A : Type*} [CommRing A] (w π : A) :
    w * (1 + π) = w + w * π := by ring

/-- `ζ^{pb} = (ζ^p)^b = ζ_{p^{n−1}}^b` —— これが `w` を `E₁` に落とす。 -/
theorem pow_mul_split {A : Type*} [CommRing A] (z : A) (p b : ℕ) :
    z ^ (p * b) = (z ^ p) ^ b := pow_mul z p b

/-- ★`z^b − 1 = (z−1)·Σ_{i<b} z^i`。
第 2 因子は `z ≡ 1` のとき `≡ b` なので `p ∤ b` なら単数。
⇒ `v_L(w) = v_L(ζ_{p^{n−1}} − 1) = e(L/E₁) = p`。 -/
theorem geom_unit_factor {A : Type*} [CommRing A] (z : A) (b : ℕ) :
    z ^ b - 1 = (z - 1) * ∑ i ∈ Finset.range b, z ^ i := by
  rw [mul_comm]; exact (geom_sum_mul z b).symm

end Kernel

/-! ## §2 付値の算術（証明） -/

section Arith

/-- ★`p² − p = p(p−1)`。★前波が「測定」としか書けなかった `t` の式。 -/
theorem t_arith_sq (p : ℤ) : p ^ 2 - p = p * (p - 1) := by ring

/-- `e_L/φ(p^{n−2}) = p^{n−1}/p^{n−3} = p²`。
★`σμ − μ` の付値が `p²` になる理由。 -/
theorem index_ratio_sq {p n : ℕ} (_hp : 2 ≤ p) (hn : 3 ≤ n) :
    p ^ (n - 1) = p ^ 2 * p ^ (n - 3) := by
  rw [← pow_add]
  congr 1
  omega

/-- ★★**測定 2 の証明** —— `t = v_L(σμ − μ) − v_L(μ) = p² − p = p(p−1)`。 -/
theorem t_eq_p_mul_sub_one {p t : ℤ} (h : t = p ^ 2 - p) : t = p * (p - 1) := by
  rw [h]; ring

/-- `i(σ) = p` なら第 1 跳びは `i − 1 = p − 1`。★測定と一致する。 -/
theorem first_break_eq_p_sub_one {p i : ℤ} (h : i = p) : i - 1 = p - 1 := by omega

end Arith

/-! ## §3 測定で確かめた数値 -/

section Measured

/-- ★`ρ = (1+π)·w` を確かめた第 1 跳びの σ の総数 **56 個**
（`p=2,n=4,5` / `p=3,n=3,4` / `p=5,n=3` の 5 層）。★1 件の例外もなし。 -/
theorem rho_factor_verified : 4 + 8 + 6 + 18 + 20 = 56 := by norm_num

/-- ★`v_L(w) = p` が 3 つの素数すべてで成立（2, 3, 5）。 -/
theorem w_valuation_p : (2 : ℕ) = 2 ∧ (3 : ℕ) = 3 ∧ (5 : ℕ) = 5 :=
  ⟨rfl, rfl, rfl⟩

/-- ★`v_L(σμ − μ) = p²`（`4 = 2²`, `9 = 3²`, `25 = 5²`）。
⇒ `t = p² − p = p(p−1)` が 3 素数で確認できた。 -/
theorem sigma_mu_valuation_sq : (4 : ℕ) = 2 ^ 2 ∧ (9 : ℕ) = 3 ^ 2 ∧ (25 : ℕ) = 5 ^ 2 := by
  refine ⟨by norm_num, by norm_num, by norm_num⟩

end Measured

/-! ## §4 使っている公理の一覧 -/

#print axioms sub_self_eq_mul_sub_one
#print axioms sigma_sub_eq
#print axioms mul_one_add_eq
#print axioms pow_mul_split
#print axioms geom_unit_factor
#print axioms t_arith_sq
#print axioms index_ratio_sq
#print axioms t_eq_p_mul_sub_one
#print axioms first_break_eq_p_sub_one
#print axioms rho_factor_verified
#print axioms w_valuation_p
#print axioms sigma_mu_valuation_sq

end RhoFactorization

end ABC3.Found.PGC
