/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.ProL.DivByP

/-!
# `u^{1-l} = v_l` を全素数で同時に解く（鎖 `prop56` の `p56-limit`）

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.107。

原文 (FrdI p.107):
> as desired, it suffices to show, for each p ∈Primes, the existence of a up ∈O×(A)[p]

## ★★何のためか

`Proposition 5.6` の証明の山は「`u · φ_l · u⁻¹ = φ′_l` をすべての素数 `l` で
同時に満たす `u ∈ 𝒪^×(A)` の存在」である。Frobenius-normalized 型だから
`u · φ_l · u⁻¹ = u^{1-l} · φ_l` であり、`φ′_l = v_l · φ_l` と書けば
（`exists_otimes_of_frobType_same_deg`）、求めるものは

  `u^{1-l} = v_l`（すべての素数 `l`）

という**連立方程式**にほかならない。本ファイルはそれを可換副有限群の言葉だけで解く。

## ★★★原典との差（逸脱の記録、2026-08-20）

原典は `E_{M_p} = M_p \ E` という剰余単系を作り、そこでの `p≈` で議論して
最後に「`u_p` の**無限積**を取る」と書く。★我々はそれを

* `M ≅ ∏_l M[l]`（`decompEquiv`、既済）で**成分ごとの方程式**に分け、
* 各 `M[q]` で `x ↦ x^{q-1}` が全単射（`pow_pred_bijective`、既済）だから
  `q` 成分の解 `u_q` を取り、
* `u := decompEquiv.symm (fun q => u_q)`

と読み替えた。★★**「無限積」は積分解の逆像そのもの**であって、
コンパクト性による逆極限の構成は要らない。
★これは原典の意図（`M_p` を潰して `p` 成分だけを見る）に忠実であり、
`Proposition 5.6` を消費する側（`Corollary 5.7`）に影響しない。

## ★連立が解ける理由（両立条件）

`φ, φ′` がともに単系準同型であることから `v` は

  `v_{l₁}^{1-l₂} = v_{l₂}^{1-l₁}`  …… (†)

を満たす（`φ′_{l₁}φ′_{l₂} = φ′_{l₂}φ′_{l₁}` を消去したもの）。
★(†) を `c_l := v_l⁻¹` で書き直すと `c_{l₁}^{l₂-1} = c_{l₂}^{l₁-1}` である。
これがちょうど「`q ≠ l` の成分でも辻褄が合う」ための条件になっている。
-/

namespace ABC3.Found.ProL

universe u

variable {M : Type u} [CommGroup M] [TopologicalSpace M] [IsTopologicalGroup M]
  [CompactSpace M] [TotallyDisconnectedSpace M]

/-- ★`q` 成分での `x ↦ x^{q-1}` の全単射性（`pow_pred_bijective` を `lPart` に当てたもの）。 -/
theorem pow_pred_bijective_lPart (q : Nat.Primes) :
    Function.Bijective (fun x : ↥(lPart M (q : ℕ)) => x ^ ((q : ℕ) - 1)) :=
  pow_pred_bijective q.2 (lPart_isProL (M := M) (q : ℕ))

/-- ★★★★★**連立 `u^{l-1} = c_l` は両立条件 (†) のもとで解ける**。

★可換副有限群 `M` と族 `c : Primes → M` が
`c l₁ ^ (l₂ - 1) = c l₂ ^ (l₁ - 1)` をすべての素数対で満たすなら、
`u ^ (l - 1) = c l` をすべての素数 `l` で満たす `u ∈ M` がある。

★★これが `Proposition 5.6` の「`u` の存在」（原文 p.107 の `u_p` の無限積）である。 -/
theorem exists_pow_pred_eq (c : Nat.Primes → M)
    (hc : ∀ l₁ l₂ : Nat.Primes, (c l₁) ^ ((l₂ : ℕ) - 1) = (c l₂) ^ ((l₁ : ℕ) - 1)) :
    ∃ u : M, ∀ l : Nat.Primes, u ^ ((l : ℕ) - 1) = c l := by
  classical
  -- ★1. 各成分で方程式を解く
  have hsol : ∀ q : Nat.Primes, ∃ y : ↥(lPart M (q : ℕ)),
      y ^ ((q : ℕ) - 1) = (decompEquiv M (c q)) q := fun q =>
    (pow_pred_bijective_lPart (M := M) q).2 _
  choose y hy using hsol
  refine ⟨(decompEquiv M).symm y, fun l => ?_⟩
  -- ★2. 成分ごとに確かめる
  refine (decompEquiv M).injective ?_
  funext q
  have hu : (decompEquiv M ((decompEquiv M).symm y)) = y :=
    (decompEquiv M).apply_symm_apply y
  rw [map_pow, hu]
  by_cases hql : q = l
  · subst hql
    exact hy q
  · -- ★3. `q ≠ l`：両立条件 (†) を `q` 成分で読む
    -- `(y q ^ (l-1)) ^ (q-1) = ((decompEquiv M (c l)) q) ^ (q-1)` を出して
    -- `x ↦ x^{q-1}` の単射性で落とす。
    refine (pow_pred_bijective_lPart (M := M) q).1 ?_
    show (y q ^ ((l : ℕ) - 1)) ^ ((q : ℕ) - 1)
      = ((decompEquiv M (c l)) q) ^ ((q : ℕ) - 1)
    have h1 : (y q ^ ((l : ℕ) - 1)) ^ ((q : ℕ) - 1)
        = (y q ^ ((q : ℕ) - 1)) ^ ((l : ℕ) - 1) := by
      rw [← pow_mul, ← pow_mul, Nat.mul_comm]
    have h2 : ((decompEquiv M (c q)) q) ^ ((l : ℕ) - 1)
        = ((decompEquiv M (c l)) q) ^ ((q : ℕ) - 1) := by
      have := congrArg (fun z : M => (decompEquiv M z) q) (hc q l)
      simpa using this
    rw [h1, hy q, h2]

/-- ★★★★**`u^{1-l} = v_l` の形**（原文の書き方に合わせたもの）。

★`u * (u ^ l)⁻¹ = v l` は `u^{1-l} = v_l` のことである。 -/
theorem exists_forall_pow_one_sub (v : Nat.Primes → M)
    (hv : ∀ l₁ l₂ : Nat.Primes,
      (v l₁) * ((v l₁) ^ (l₂ : ℕ))⁻¹ = (v l₂) * ((v l₂) ^ (l₁ : ℕ))⁻¹) :
    ∃ u : M, ∀ l : Nat.Primes, u * (u ^ (l : ℕ))⁻¹ = v l := by
  -- ★`x * (x^n)⁻¹ = (x^{n-1})⁻¹`（`n ≥ 1`）に書き換えてから `exists_pow_pred_eq`。
  have key : ∀ (x : M) (n : ℕ), 1 ≤ n → x * (x ^ n)⁻¹ = (x ^ (n - 1))⁻¹ := by
    intro x n hn
    have hsplit : x ^ n = x ^ (n - 1) * x := by
      rw [← pow_succ]; congr 1; omega
    rw [hsplit, mul_inv_rev]
    group
  have hc : ∀ l₁ l₂ : Nat.Primes,
      ((v l₁)⁻¹) ^ ((l₂ : ℕ) - 1) = ((v l₂)⁻¹) ^ ((l₁ : ℕ) - 1) := by
    intro l₁ l₂
    have e := hv l₁ l₂
    rw [key _ _ l₂.2.one_lt.le, key _ _ l₁.2.one_lt.le] at e
    rw [inv_pow, inv_pow, e]
  obtain ⟨u, hu⟩ := exists_pow_pred_eq (fun l => (v l)⁻¹) hc
  refine ⟨u, fun l => ?_⟩
  rw [key _ _ l.2.one_lt.le, hu l, inv_inv]

/-! ### ★出典の紐付け -/

/-- ★★locator —— `Proposition 5.6` の「`u` の存在」（原文 p.107 の `u_p` の無限積）。 -/
def exists_forall_pow_one_sub.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 107,
    item := "Proposition 5.6 — u^{1-l} = v_l を全素数で同時に解く",
    sectionId := "frdi-prop-5-6" }

end ABC3.Found.ProL
