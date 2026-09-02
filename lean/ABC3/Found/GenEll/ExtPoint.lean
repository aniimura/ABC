/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.SSCurveExt
import ABC3.Found.GenEll.DegreeLine
import ABC3.Meta.Claim

/-!
# 第 1346 ブロック —— **降ろした点を `ℂ` の中へ運ぶ**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か——`E.alg` の中の `M` を `ℂ` の中へ

`exists_point_descent_of_stable`（第 1219）が与える中間体は
`M : IntermediateField E.fld E.alg` である。
★一方 `SSCurve.ext`（第 1343）が要るのは `IntermediateField E.fld ℂ` である。

☆`E.alg` は `E.fld` 上代数的で `ℂ` は代数閉体だから、
`ι : E.alg →ₐ[E.fld] ℂ`（`IsAlgClosed.lift`）で運べる。
★`M.map ι` が求める中間体であり、点は `equivMap` ＋ `rhPoint` で移る。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve IsDedekindDomain NumberField
open ABC3.Found.GaloisRep ABC3.Interface.GaloisRep ABC3.Meta

open scoped Classical

/-- ★★★★★★**代数閉包から `ℂ` への `E.fld`-埋め込み**（第 1346）。 -/
noncomputable def algToComplex (E : SSCurve) : E.alg →ₐ[E.fld] ℂ :=
  IsAlgClosed.lift (R := E.fld) (S := E.alg) (M := ℂ)

/-- ★★★★★★★★★★★★★★★★★★★★
**`Gal`-安定な点は `ℂ` の中の有限拡大の上で有理になる**——★**無条件**（第 1346）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆第 1219 の降下で `M : IntermediateField E.fld E.alg` を取り、
`algToComplex` で `ℂ` の中へ運ぶ。★点は `equivMap` ＋ `rhPoint` で移る。 -/
theorem exists_ext_point_of_stable (E : SSCurve) {l : ℕ} (hl : l.Prime)
    (Q : (E.W.baseChange E.alg).toAffine.Point) (hQ : addOrderOf Q = l)
    (hst : ∀ σ : E.alg ≃ₐ[E.fld] E.alg, ∃ k : ℕ, galPoint E.W σ Q = k • Q) :
    ∃ (M₀ : IntermediateField E.fld ℂ), FiniteDimensional E.fld M₀ ∧
      Module.finrank E.fld M₀ ≤ l - 1 ∧
      letI : DecidableEq ((M₀ : IntermediateField E.fld ℂ) : Type) :=
        fun a b => Classical.propDecidable (a = b)
      ∃ Q₀ : (E.W.baseChange ((M₀ : IntermediateField E.fld ℂ) : Type)).toAffine.Point,
        addOrderOf Q₀ = l := by
  classical
  obtain ⟨M, hfin, hdeg, Q', hQ', -⟩ :=
    exists_point_descent_of_stable E.W hl Q hQ hst
  haveI := hfin
  letI : DecidableEq (M : Type) := fun a b => Classical.propDecidable (a = b)
  set ι : E.alg →ₐ[E.fld] ℂ := algToComplex E with hι
  letI : DecidableEq ((M.map ι : IntermediateField E.fld ℂ) : Type) :=
    fun a b => Classical.propDecidable (a = b)
  set e : M ≃ₐ[E.fld] (M.map ι) := IntermediateField.equivMap M ι with he
  haveI hfin₀ : FiniteDimensional E.fld (M.map ι) :=
    LinearEquiv.finiteDimensional e.toLinearEquiv
  refine ⟨M.map ι, hfin₀, ?_, ?_⟩
  · rw [← LinearEquiv.finrank_eq e.toLinearEquiv]
    exact hdeg
  · -- ★点を移す
    have hcomp : (e.toAlgHom.toRingHom).comp (algebraMap E.fld M)
        = algebraMap E.fld (M.map ι) :=
    RingHom.ext fun y => e.commutes y
    have hcurve : (E.W.baseChange (M : Type)).map (e.toAlgHom.toRingHom)
        = E.W.baseChange ((M.map ι : IntermediateField E.fld ℂ) : Type) := by
      show (E.W.map (algebraMap E.fld M)).map (e.toAlgHom.toRingHom) = _
      rw [WeierstrassCurve.map_map, hcomp]
      rfl
    haveI hell1 : (E.W.baseChange (M : Type)).IsElliptic := by
      show (E.W.map (algebraMap E.fld M)).IsElliptic
      infer_instance
    haveI hell2 : ((E.W.baseChange (M : Type)).map (e.toAlgHom.toRingHom)).IsElliptic := by
      rw [hcurve]
      show (E.W.map (algebraMap E.fld (M.map ι))).IsElliptic
      infer_instance
    refine ⟨hcurve ▸ rhPoint (e.toAlgHom.toRingHom) (E.W.baseChange (M : Type)) Q', ?_⟩
    have hord := addOrderOf_congr_curve hcurve
      (rhPoint (e.toAlgHom.toRingHom) (E.W.baseChange (M : Type)) Q')
    rw [hord, addOrderOf_rhPoint]
    exact hQ'

/-! ## ★出典の紐付け(`.src`) -/

def algToComplex.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(代数閉包から ℂ への E.fld-埋め込み)",
    sectionId := "genell-lemma-3-5" }

def exists_ext_point_of_stable.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Gal-安定な点は ℂ の中の有限拡大の上で有理になる。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_ext_point_of_stable.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_point_descent_of_stable(第 1219、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_point_descent_of_stable") 1,
    .citation "[ABC3]" "SSCurve.ext(第 1343、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.SSCurve.ext") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1346）**——`IsQuotClassJ` の存在を出すための" ++
       "**運搬の段**である。☆`E.alg` の中の中間体を `ℂ` の中へ移す。") 2 ]

end ABC3.Found.GenEll
