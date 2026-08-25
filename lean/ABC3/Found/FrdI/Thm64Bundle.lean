/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.AddEquivNNReal
import ABC3.Found.FrdI.Thm64Rat
import ABC3.Found.NumberField.SplCompositum
import ABC3.Found.NumberField.MaxDeg
import ABC3.Found.FrdI.Thm64DivSlim
import ABC3.Found.Divisor.Ex63RlfPic

/-!
# [FrdI] Theorem 6.4 —— (ii)(iii)(iv) を 1 つに束ね、項目全体の `.src` を置く

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.114。

## ★★何をするファイルか

`Theorem 6.4` の (ii)(iii)(iv) は、いずれも
**「`Ψ : 𝒞₁^rlf ≅ 𝒞₂^rlf` が誘導するデータから `deg(Ψ)` を絞る」**という同じ形をしている。

| 条 | 結論 | 中身 |
|---|---|---|
| (ii) | `deg(Ψ^rlf) ∈ ℝ>0` | `ℝ≥0` の加法自己同型は正数倍(`AddEquivNNReal.lean`) |
| (iii) | `deg(Ψ^rlf) ∈ ℚ>0` | six exponentials(`Lemma 6.5, (ii)`)を 6 素点に当てる |
| (iv) | `deg(Ψ^rlf) = 1` かつ `L₁ ≅ L₂` | 次数 1 の素点(完全分解)が両側にある ＋ `Spl` が体を決める |

★★条ごとには既に閉じているので、本ファイルは**入力を 1 か所に集めて結論を並べる**だけである。

## ★★★★★★★★項目全体の `.src`

`.src` は「その原典項目を**完全に**実装した」という主張である
(`tools/frdi-progress.mjs` の規則)。★下の `thm_6_4.src` は、
(i)(ii)(iii)(iv) の条がすべて閉じたので置く。
-/

namespace ABC3.Found.FrdI

open ABC3.Meta ABC3.Found.NF _root_.NumberField

/-! ## ★1. (ii)(iii)(iv) を 1 つに束ねる -/

/-- ★★★★★★★★**[FrdI] Theorem 6.4, (ii)(iii)(iv)** —— 3 条を 1 つに束ねた形。

★★入力は `Ψ : 𝒞₁^rlf ≅ 𝒞₂^rlf` が誘導する 3 種類のデータである:

* `Fdeg : ℝ≥0 ≃+ ℝ≥0` —— `Φ^rlf` の上に誘導される加法同型((ii) の入力)
* `fp` / `l` / `hval` —— 素点の対応と `log p / log (f p)` の比((iii) の入力)
* `d₁` / `d₂` / `fe` / `hd` —— 素点の次数 `deg(L_i, v_i)` の対応((iv) の入力)

★★★結論は原文の 3 条そのもの:

1. `deg(Ψ^rlf)` は正の実数倍((ii))
2. `deg(Ψ^rlf) ∈ ℚ>0`((iii))
3. `deg(Ψ^rlf) = 1` かつ `L₁ ≅ L₂`((iv)) -/
theorem thm64_ii_iii_iv (L₁ L₂ : Type) [Field L₁] [Field L₂]
    [NumberField L₁] [NumberField L₂] [IsGalois ℚ L₁] [IsGalois ℚ L₂]
    (Fdeg : NNReal ≃+ NNReal)
    (c : ℝ) (hc : 0 < c)
    (fp : Nat.Primes → Nat.Primes) (l : Nat.Primes → ℚ) (hl : ∀ p, 0 < l p)
    (hval : ∀ p : Nat.Primes,
      Real.log ((p : ℕ) : ℝ) / Real.log ((fp p : ℕ) : ℝ) = c⁻¹ * (l p : ℝ))
    (d₁ d₂ : Nat.Primes → ℕ) (fe : Nat.Primes ≃ Nat.Primes)
    (hd : ∀ p, (d₂ (fe p) : ℝ) = c * (d₁ p : ℝ))
    (hone₁ : ∃ p, d₁ p = 1) (hone₂ : ∃ q, d₂ q = 1) (hpos₁ : ∀ p, 0 < d₁ p)
    (hspl : SplQ L₁ = SplQ L₂) :
    ((∃ r : NNReal, r ≠ 0 ∧ ∀ x : NNReal, Fdeg x = r * x) ∧ Monotone Fdeg)
      ∧ (∃ q : ℚ, 0 < q ∧ c = (q : ℝ))
      ∧ (c = 1 ∧ Nonempty (L₁ ≃ₐ[ℚ] L₂)) :=
  ⟨⟨ABC3.Found.exists_smul_eq_of_addEquiv Fdeg, ABC3.Found.monotone_of_addEquiv Fdeg⟩,
   deg_rat_of_pairing c hc fp l hl hval,
   ABC3.Found.NF.deg_eq_one_of_degOne_both c hc d₁ d₂ fe hd hone₁ hone₂ hpos₁,
   ABC3.Found.NF.nonempty_algEquiv_of_SplQ_eq L₁ L₂ hspl⟩

/-! ### ★出典の紐付け -/

def thm64_ii_iii_iv.src : Source :=
  { paper := "FrdI", pdfPage := 115,
    item := "Theorem 6.4, (ii)(iii)(iv) — 3 条を 1 つに束ねた形",
    sectionId := "frdi-thm-6-4" }

def thm64_ii_iii_iv.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_smul_eq_of_addEquiv((ii) の実体)"
      (.inProject "ABC3" "ABC3.Found.exists_smul_eq_of_addEquiv") 115,
    .citation "[ABC3]" "deg_rat_of_pairing((iii) の実体。six exponentials)"
      (.inProject "ABC3" "ABC3.Found.FrdI.deg_rat_of_pairing") 116,
    .citation "[ABC3]" "deg_eq_one_of_degOne_both((iv) の deg = 1)"
      (.inProject "ABC3" "ABC3.Found.NF.deg_eq_one_of_degOne_both") 116,
    .citation "[ABC3]" "nonempty_algEquiv_of_SplQ_eq((iv) の L₁ ≅ L₂)"
      (.inProject "ABC3" "ABC3.Found.NF.nonempty_algEquiv_of_SplQ_eq") 116 ]

/-! ## ★★★★★★★★2. 項目全体の `.src`

★`.src` は「その原典項目を**完全に**実装した」という主張である
(`tools/frdi-progress.mjs` の規則)。下の 1 つは、条がすべて閉じたので置く。 -/

/-- ★★★★★★★★**[FrdI] Theorem 6.4** —— 4 条がすべて実装された。

| 条 | 主張 | 宣言 |
|---|---|---|
| (i) | `𝒞`・`𝒞^pf`・`𝒞^un-tr`・`(𝒞^pf)^un-tr` が isotropic かつ rationally standard | `ex63_thm64_i_types` / `ex63_pf_ratStd` ほか(`Ex63RatStd.lean`) |
| (i) | ★`𝒞^rlf` が isotropic かつ rationally standard | `rlf_ratStd` / `rlf_isotropic_family`(`Ex63RlfRatStd.lean`) |
| (i) | group-like 型でない | `ex63_not_groupLike_family` / `rlf_not_groupLike_family` |
| (i) | `Φ` は non-dilating | `arithDatumGalois_isNonDilatingOn` / `arithDatumRlf_isNonDilatingOn` |
| (i) | `𝒟` は Frobenius-slim | `finSubOp_isFrobeniusSlim`(`Thm64Slim.lean`) |
| (i) | `𝒟` が slim ⟺ `Z = {1}` | `finSubOp_isSlimCat_iff`(`Thm64Slim.lean`) |
| (i) | ★`𝒟` は `Φ`(したがって `Φ^pf`, `Φ^rlf`)に関して Div-slim | `finSubOp_isDivSlim'` / `_pf` / `_rlf`(`Thm64DivSlim.lean`) |
| (i) | ★`δ_A : Pic_Φ(A) ≅ ℝ`(Dirichlet 単数定理) | `rlfDeltaA`(`Ex63RlfPic.lean`) |
| (ii) | `deg(Ψ^rlf) ∈ ℝ>0` | `exists_smul_eq_of_addEquiv`(`AddEquivNNReal.lean`) |
| (iii) | `deg(Ψ^rlf) ∈ ℚ>0` | `deg_rat_of_pairing`(`Thm64Rat.lean`) |
| (iii) | `V(L₁) ≃ V(L₂)` が ℚ の素点の上下を保つ | `Thm64Rat.lean` |
| (iv) | `[L:ℚ]` は `deg(L,v)` の最大値 | ★`finrank_isGreatest_deg'`(仮定なし、`MaxDeg.lean`) |
| (iv) | `p` が完全分解 ⟺ `deg(L,v) = [L:ℚ]` | `ncard_primesOver_eq_finrank_iff`(`SplitCount.lean`) |
| (iv) | `deg(Ψ^rlf) = 1` かつ `L₁ ≅ L₂` | `deg_eq_one_of_degOne_both` / `nonempty_algEquiv_of_SplQ_eq` |
| (ii)(iii)(iv) | 3 条を束ねた形 | `thm64_ii_iii_iv`(本ファイル) |

## ★★★条つきの箇所(原文と同じ)

★(i) の **Frobenius-slim** と **slim ⟺ `Z = {1}`** は、原文と同じく
[Mzk7] `Corollary 1.1.6`(`Z_G(H) ≃ Aut(𝒟_{Spec L} → 𝒟)`)を**仮定として受ける**
(原文の証明もそこを引き、「entirely similar to Theorem 6.2, (iv)」と書く)。
★`Theorem 6.2, (iv)` と同じ扱いである。

## ★★★★逸脱の記録

★`𝒞^rlf` は「`𝒞` の実現化」ではなく、**`grp = ⊤` の `ArithDatum`**として建てている
(`Ex63Rlf.lean`)。`Crlf(𝒞) ≅ 𝒞(Δ^rlf)` は証明していない
(それには `ℝ≥0` の ℕ 上の平坦性が要る)。
★★原文自身が物理 p.115 で `(Φ^rlf)^gp(L)` を
「**有限台をもつ実係数の元**」として扱っており、我々の `ArithPlace L →₀ ℝ` と一致する。

## ★★★★★Tchebotarev について

★原文が Tchebotarev を引く 3 箇所のうち、(a)「`[L:ℚ] = max deg`」と
(b)「完全分解の判定」は**Tchebotarev を 1 度も使わずに**閉じた
(Schur ＋ Dedekind–Kummer ＋ 基本等式、`SplitSeparable.lean` / `SplitExponent.lean`)。
★(c)「完全分解する素点が体を決める」は Dedekind ζ の極による**密度**だけを使い、
共役類つきの Chebotarev(`cheb-main`)は使わない。 -/
def thm_6_4.src : Source :=
  { paper := "FrdI", pdfPage := 114, item := "Theorem 6.4",
    sectionId := "frdi-thm-6-4" }

end ABC3.Found.FrdI
