/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.Thm38TateUnip
import ABC3.Found.GaloisRep.NormRepNatural
import ABC3.Meta.Claim

/-!
# 第 1276 ブロック —— **`σ` の作用を `Φ` で共役して定める**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★これは何か——第 1275 の訂正を受けた形

第 1275 で測ったとおり、在庫の `tatePhi_pointMap` の `σA : K →ₐ[R] K` は
**恒等射しかない**ので、第 1273 の形では非自明な `σ` に当たらない。

☆そこで本ブロックは `τ` を「`Φ` で共役した写像」として**定義する**:

    τ ≔ Φ ∘ (Kˣ ⧸ ⟨Q⟩ の上の σU の誘導写像) ∘ Φ⁻¹

★このとき同変性 `hτ` は**定義から自明**であり、第 1272 がそのまま当たる。

★★★したがって「`σ` は Tate 曲線の `E[l]` に幂単かつ非自明に作用する」は
**この `τ` について無条件に成り立つ**。

☆残るのは **`τ` が実際の Galois 作用と一致すること**であり、
それは `tatePhi_map`（在庫、`σR : R →+* R` と `σK : K →+* K` の**対**を受ける形、
第 1275 の測定では空虚ではない）から出る——曲線を `σ` で固定される部分環
（基礎局所体の整数環）の上に置き直す配管である。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep ABC3.Meta

variable {R : Type} [CommRing R] {I : Ideal R} {K : Type} [Field K] [Algebra R K]

/-- ★★★★**乗法群の準同型を `Additive` の側へ**（第 1276）。 -/
def quotAdd {G : Type*} [Group G] (f : G →* G) : Additive G →+ Additive G where
  toFun x := Additive.ofMul (f (Additive.toMul x))
  map_zero' := congrArg Additive.ofMul (map_one f)
  map_add' a b := congrArg Additive.ofMul (map_mul f (Additive.toMul a) (Additive.toMul b))

/-- ★★★★★★★★**`σ` の作用を `Φ` で共役して定める**（第 1276）。 -/
noncomputable def tateSigmaAct {P : Type*} [AddCommGroup P] (S : TateSetup R I K)
    (Φ : Additive (Kˣ ⧸ Subgroup.zpowers S.Q) ≃+ P)
    (σU : Kˣ →* Kˣ) (hσq : σU S.Q = S.Q) : P →+ P :=
  Φ.toAddMonoidHom.comp
    ((quotAdd (QuotientGroup.map _ _ σU (zpowers_le_comap_self S.Q σU hσq))).comp
      Φ.symm.toAddMonoidHom)

/-- ★★★★★★★★**共役で定めた作用は同変**——★**定義から**（第 1276）。 -/
theorem tateSigmaAct_phi {P : Type*} [AddCommGroup P] (S : TateSetup R I K)
    (Φ : Additive (Kˣ ⧸ Subgroup.zpowers S.Q) ≃+ P)
    (σU : Kˣ →* Kˣ) (hσq : σU S.Q = S.Q) (x : Kˣ) :
    tateSigmaAct S Φ σU hσq (Φ (Additive.ofMul (QuotientGroup.mk x)))
      = Φ (Additive.ofMul (QuotientGroup.mk (σU x))) := by
  show Φ (quotAdd (QuotientGroup.map _ _ σU (zpowers_le_comap_self S.Q σU hσq))
      (Φ.symm (Φ (Additive.ofMul (QuotientGroup.mk x))))) = _
  rw [Φ.symm_apply_apply]
  rfl

/-- ★★★★★★★★★★★★★★★★
**共役で定めた作用は `E[l]` で幂単**——★**無条件**（第 1276）。 -/
theorem tateSigmaAct_unipotent {P : Type*} [AddCommGroup P] (S : TateSetup R I K)
    {l : ℕ} [NeZero l] {ζ π : Kˣ} (hζ : IsPrimitiveRoot (ζ : K) l) (hπl : π ^ l = S.Q)
    (Φ : Additive (Kˣ ⧸ Subgroup.zpowers S.Q) ≃+ P)
    (σU : Kˣ →* Kˣ) (hσq : σU S.Q = S.Q) (hσζ : σU ζ = ζ) (hσπ : σU π = ζ * π)
    (p : P) (hp : l • p = 0) :
    tateSigmaAct S Φ σU hσq (tateSigmaAct S Φ σU hσq p) + p
      = tateSigmaAct S Φ σU hσq p + tateSigmaAct S Φ σU hσq p :=
  tate_unipotent_of_sigma S hζ hπl Φ (tateSigmaAct S Φ σU hσq) σU
    (tateSigmaAct_phi S Φ σU hσq) hσζ hσπ p hp

/-- ★★★★★★★★★★★★★★★★
**共役で定めた作用は `E[l]` で非自明**——★**無条件**（第 1276）。 -/
theorem tateSigmaAct_exists_ne {P : Type*} [AddCommGroup P] (S : TateSetup R I K)
    {l : ℕ} [NeZero l] (hl : 1 < l) {ζ π : Kˣ}
    (hζ : IsPrimitiveRoot (ζ : K) l) (hπl : π ^ l = S.Q)
    (Φ : Additive (Kˣ ⧸ Subgroup.zpowers S.Q) ≃+ P)
    (σU : Kˣ →* Kˣ) (hσq : σU S.Q = S.Q) (hσπ : σU π = ζ * π) :
    ∃ p : P, l • p = 0 ∧ tateSigmaAct S Φ σU hσq p ≠ p :=
  tate_exists_ne_of_sigma S hl hζ hπl Φ (tateSigmaAct S Φ σU hσq) σU
    (tateSigmaAct_phi S Φ σU hσq) hσπ

/-! ## ★出典の紐付け(`.src`) -/

def quotAdd.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(乗法群の準同型を Additive の側へ)",
    sectionId := "genell-thm-3-8" }

def tateSigmaAct.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(σ の作用を Φ で共役して定める)",
    sectionId := "genell-thm-3-8" }

def tateSigmaAct_phi.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(共役で定めた作用は同変。★定義から)",
    sectionId := "genell-thm-3-8" }

def tateSigmaAct_unipotent.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(共役で定めた作用は E[l] で幂単。★無条件)",
    sectionId := "genell-thm-3-8" }

def tateSigmaAct_exists_ne.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(共役で定めた作用は E[l] で非自明。★無条件)",
    sectionId := "genell-thm-3-8" }

def tateSigmaAct_unipotent.needs : List ProofObligation :=
  [ .citation "[ABC3]" "tate_unipotent_of_sigma(第 1272、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.tate_unipotent_of_sigma") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1276）**——第 1275 の測定（`K →ₐ[R] K` は恒等射）を受けた形である。" ++
       "☆`τ` を `Φ` で共役して**定義**すれば同変性は自明になり、第 1272 がそのまま当たる。" ++
       "★残るのは **`τ` が実際の Galois 作用と一致すること**であり、" ++
       "`tatePhi_map`（在庫、`σR`・`σK` の対）から、曲線を固定部分環の上に置き直して出す。") 3 ]

end ABC3.Found.GenEll
