/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.Thm38SigmaExists
import Mathlib.RingTheory.RootsOfUnity.PrimitiveRoots
import ABC3.Meta.Claim

/-!
# 第 1281 ブロック —— **`σ(π) = ζπ` で `ζ` は原始 `l` 乗根**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★これは何か——`α` の形が出る段

第 1280 で `π` を動かす `σ` を取った。★本ブロックはその `σ` が
**ちょうど `α = (1 1 / 0 1)` の形で作用する**ことを言う:

> `σ(π)^l = σ(Q) = Q = π^l` なので `ζ ≔ σ(π)/π` は `l` 乗根。
> `σ(π) ≠ π` なので `ζ ≠ 1`、`l` は素数なので `ζ` は**原始 `l` 乗根**。

★★★これで `tateSigmaAct_unipotent`・`tateSigmaAct_exists_ne`（第 1276）が
要求する `hσπ : σU π = ζ * π` と `hζ : IsPrimitiveRoot ζ l` が**同時に**出る。
-/

namespace ABC3.Found.GenEll

open ABC3.Meta

/-- ★★★★★★★★★★★★★★★★★★★★
**`σ` が `π` を動かせば `σ(π) = ζπ` で `ζ` は原始 `l` 乗根**——★**無条件**（第 1281）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`ζ ≔ σ(π)/π` の `l` 乗は `σ(Q)/Q = 1`。★`ζ ≠ 1` と `l` 素数から位数はちょうど `l`。 -/
theorem sigma_pi_eq_zeta_mul {K₀ M : Type} [Field K₀] [Field M] [Algebra K₀ M]
    {l : ℕ} (hl : l.Prime) (Q₀ : K₀ˣ) (π : Mˣ)
    (hπ : π ^ l = Units.map (algebraMap K₀ M : K₀ →* M) Q₀)
    (σ : M ≃ₐ[K₀] M) (hne : σ (π : M) ≠ (π : M)) :
    ∃ ζ : Mˣ, IsPrimitiveRoot ((ζ : M)) l ∧
      Units.map σ.toAlgHom.toRingHom.toMonoidHom π = ζ * π := by
  set σM : M →* M := σ.toAlgHom.toRingHom.toMonoidHom with hσM
  set ζ : Mˣ := Units.map σM π * π⁻¹ with hζdef
  have hfix : Units.map σM (Units.map (algebraMap K₀ M : K₀ →* M) Q₀)
      = Units.map (algebraMap K₀ M : K₀ →* M) Q₀ := by
    refine Units.ext ?_
    show σ (algebraMap K₀ M ((Q₀ : K₀))) = algebraMap K₀ M ((Q₀ : K₀))
    exact σ.commutes _
  have hζl : ζ ^ l = 1 := by
    have h1 : (Units.map σM π) ^ l = Units.map σM (π ^ l) := (map_pow _ _ _).symm
    have h2 : (π⁻¹ : Mˣ) ^ l = (π ^ l)⁻¹ := inv_pow π l
    rw [hζdef, mul_pow, h1, h2, hπ, hfix, mul_inv_cancel]
  have hζne : ζ ≠ 1 := by
    intro h
    refine hne ?_
    have hmap : Units.map σM π = π := by
      have h3 := congrArg (fun z : Mˣ => z * π) h
      simpa [hζdef, inv_mul_cancel_right] using h3
    exact congrArg (fun u : Mˣ => (u : M)) hmap
  have hord : orderOf ζ = l := by
    have h1 : orderOf ζ ∣ l := orderOf_dvd_of_pow_eq_one hζl
    rcases (Nat.Prime.eq_one_or_self_of_dvd hl _ h1) with h | h
    · exact absurd (orderOf_eq_one_iff.mp h) hζne
    · exact h
  refine ⟨ζ, ⟨?_, ?_⟩, ?_⟩
  · have h4 : ((ζ ^ l : Mˣ) : M) = 1 := by rw [hζl]; rfl
    simpa using h4
  · intro m hm
    have hu : ζ ^ m = 1 := by
      refine Units.ext ?_
      push_cast
      simpa using hm
    have hd := orderOf_dvd_of_pow_eq_one hu
    rwa [hord] at hd
  · rw [hζdef, inv_mul_cancel_right]

/-! ## ★出典の紐付け(`.src`) -/

def sigma_pi_eq_zeta_mul.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(σ が π を動かせば σ(π) = ζπ で ζ は原始 l 乗根。★無条件)",
    sectionId := "genell-thm-3-8" }

def sigma_pi_eq_zeta_mul.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_algEquiv_move_pi(第 1280、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_algEquiv_move_pi") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1281）**——`α = (1 1 / 0 1)` の形が出る段である。" ++
       "☆`tateSigmaAct_unipotent`・`tateSigmaAct_exists_ne`（第 1276）が要求する" ++
       "`hσπ : σU π = ζ * π` と `hζ : IsPrimitiveRoot ζ l` が**同時に**出る。") 2 ]

end ABC3.Found.GenEll
