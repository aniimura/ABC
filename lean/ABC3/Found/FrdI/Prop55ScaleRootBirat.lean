/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop55ScaleRootCoa

/-!
# [FrdI] Proposition 5.5, (ii) —— `Σ_k` の birat 版は底と次数を保つ

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.105。

原文 (FrdI p.105):
> if C is of standard (respectively, rationally standard) type, then so are Cun-tr, Crlf.

## ★★何のために要るか

`Theorem 6.4, (i)` の 5 圏のうち `𝒞^pf` を止めているのは
**`𝒞^pf` が birationally Frobenius-normalized 型**という 1 条である
(`Skeleton/FrdI/Prop55PfBiratFn.lean` に主張と手順書がある)。

★その手順書の第 3 段が「`Σ_k` の birat 版(`scaleRootBirat`)の
`Base`・`degFr` の保存」であり、**本ファイルがそれを埋める**。

## ★★中身は代表元に降ろすだけ

`scaleRootBirat` は `psiBiratCor`(`Corollary 4.10`)で作ってあるので、
`biratPsiMap_mk`(代表元での値)で `HomBirat.mk Z φ` に降ろせば、
あとは **`𝒞^pf` の側の保存**(在庫の `rootBase_scaleRootHom` /
`rootDeg_scaleRootHom`)がそのまま効く。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section ScaleRootBirat

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}
  [IsConnected D]

/-- ★★★★**`Σ_k` の birat 版は Frobenius 次数を保つ**。 -/
theorem biratDeg_scaleRootBirat (k : ℕ+) (Gpf : Frobenioid (pfRootPre P F))
    {X Y : BiratCat (pfRootPre P F) Gpf} (f : X ⟶ Y) :
    biratDeg ((scaleRootBirat (F := F) k Gpf).map f) = biratDeg f := by
  obtain ⟨Z, φ, hZ⟩ := HomBirat.exists_rep f
  subst hZ
  show biratDeg (biratPsiMap Gpf Gpf (scaleRootFunctor P F k)
      (fun {_ _} g hg => coaPreProp_scaleRootHom k g hg) X Y (HomBirat.mk Z φ)) = _
  rw [biratPsiMap_mk, biratDeg_mk, biratDeg_mk]
  exact rootDeg_scaleRootHom k φ

/-- ★★★★★**`Σ_k` の birat 版は底を保つ**。

★`biratBase` は `sliceBaseOf a _ φ = inv (Base a) ≫ Base φ` なので、
`Base` が両方とも保たれる(`rootBase_scaleRootHom`)ことから出る。 -/
theorem biratBase_scaleRootBirat (k : ℕ+) (Gpf : Frobenioid (pfRootPre P F))
    {X Y : BiratCat (pfRootPre P F) Gpf} (f : X ⟶ Y) :
    biratBase ((scaleRootBirat (F := F) k Gpf).map f) = biratBase f := by
  obtain ⟨Z, φ, hZ⟩ := HomBirat.exists_rep f
  subst hZ
  show biratBase (biratPsiMap Gpf Gpf (scaleRootFunctor P F k)
      (fun {_ _} g hg => coaPreProp_scaleRootHom k g hg) X Y (HomBirat.mk Z φ)) = _
  rw [biratPsiMap_mk, biratBase_mk, biratBase_mk]
  have hb1 : (pfRootPre P F).Base
      ((scaleRootFunctor P F k).map Z.unop.hom.hom)
      = (pfRootPre P F).Base Z.unop.hom.hom :=
    rootBase_scaleRootHom k Z.unop.hom.hom
  have hb2 : (pfRootPre P F).Base ((scaleRootFunctor P F k).map φ)
      = (pfRootPre P F).Base φ := rootBase_scaleRootHom k φ
  have hidx : (Over.hom (Opposite.unop
      ((idxBiratPsi Gpf Gpf (scaleRootFunctor P F k)
        (fun {_ _} g hg => coaPreProp_scaleRootHom k g hg) X).obj Z))).hom
      = (scaleRootFunctor P F k).map Z.unop.hom.hom := rfl
  have hbase := (congrArg (pfRootPre P F).Base hidx).trans hb1
  simp only [sliceBaseOf, hbase]
  congr 1

end ScaleRootBirat

/-! ### ★出典の紐付け -/

def biratDeg_scaleRootBirat.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — Σ_k の birat 版は Frobenius 次数を保つ",
    sectionId := "frdi-prop-5-5" }

def biratBase_scaleRootBirat.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — Σ_k の birat 版は底を保つ",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
