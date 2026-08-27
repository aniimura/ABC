/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ConductorDescent

/-!
# イデアル層の引き戻しを**アフィンで計算する**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

## ★★★★★★★★★★`mathlib-gap` の `ideal-comap-on-affine-opens` の**アフィンの場合**

台帳（`ResearchPaper/mathlib-gap.json`）はこう書いていた:

> 一般の射 `f : X ⟶ Y` とイデアル層 `I` に対し、アフィン開 `U ⊆ f⁻¹(V)` の上で
> `(I.comap f).ideal U = (I.ideal V).map (f.appLE V U)` と書く補題。★mathlib に**無い**。

★★本ファイルは **`X`・`Y` がともにアフィンの場合**を閉じる。

    (ker (Spec R → Spec (R⧸J)) を Spec S へ引き戻したもの).ideal ⊤ = J·S

## ★★★★★★組み立ての 3 本

| 補題 | 中身 | 道具 |
|---|---|---|
| `ker_specMap_ideal_top` | `ker(Spec.map ψ).ideal ⊤ = comap ΓSpecIso (ker ψ)` | `Scheme.Hom.ker_apply`・`ΓSpecIso_naturality` |
| `ker_algebraMap_tensor` | `ker(S → S ⊗_R R⧸J) = J·S` | `quotIdealMapEquivTensorQuot` |
| `comap_ker_quotient` | 引き戻しは `Spec(S ⊗_R R⧸J)` の核 | `ker_fst_of_isClosedImmersion`・`pullbackSpecIso_inv_fst′` |

★★★**部品はすべて mathlib にあった**——無かったのは繋ぎ方だけである。

## ★★★★これで何が開くか

★`comap_mul`（`(I*J).comap f = (I.comap f)*(J.comap f)`）は
本補題と `Ideal.map_mul` から出る。mathlib には `comap_mono` はあるが
`comap_mul` は**無い**（2026-08-27 実測）。

★★その先に `Remark 1.4.1` の段 3b（チャートごとの `m` を有限被覆で揃える）がある。

## ★残っている段（明示）

★**一般の（アフィンでない）`X`・`Y` の場合**。
道筋は決まっている——`comap_comp` で開埋め込みに沿って制限し、
`ideal_comap_of_isOpenImmersion`（mathlib に**ある**）で `.ideal U` に降ろし、
本ファイルのアフィンの場合を当てる。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Limits TensorProduct

/-! ## ★`Spec.map` の核をアフィンで読む -/

/-- ★★★**`ker (Spec.map ψ)` の `⊤` での値は `ker ψ` の引き戻し**。

★`Scheme.Hom.ker_apply`（`[QuasiCompact]` が要る——`Spec.map` は常にアフィン射）と
`ΓSpecIso_naturality` を繋ぐだけ。

★★`Γ(Spec A, ⊤)` は `A` と同型だが**同じ型ではない**ので、
`ΓSpecIso` の `comap` が入る。 -/
theorem ker_specMap_ideal_top {A B : CommRingCat} (ψ : A ⟶ B) :
    (Scheme.Hom.ker (Spec.map ψ)).ideal ⟨⊤, AlgebraicGeometry.isAffineOpen_top (Spec A)⟩
      = Ideal.comap (Scheme.ΓSpecIso A).hom.hom (RingHom.ker ψ.hom) := by
  have happ : Scheme.Hom.appTop (Spec.map ψ)
      = ((Scheme.ΓSpecIso A).hom ≫ ψ) ≫ (Scheme.ΓSpecIso B).inv :=
    (Iso.eq_comp_inv _).2 (Scheme.ΓSpecIso_naturality ψ)
  have hinj : Function.Injective (Scheme.ΓSpecIso B).inv.hom :=
    ConcreteCategory.injective_of_mono_of_preservesPullback (Scheme.ΓSpecIso B).inv
  rw [Scheme.Hom.ker_apply]
  show RingHom.ker (CommRingCat.Hom.hom (Scheme.Hom.appTop (Spec.map ψ))) = _
  rw [happ]
  ext x
  simp only [CommRingCat.hom_comp, RingHom.mem_ker, RingHom.coe_comp, Function.comp_apply,
    Ideal.mem_comap]
  constructor
  · intro hx
    exact hinj (by rw [hx, map_zero])
  · intro hx
    rw [hx, map_zero]

/-! ## ★★環の側 —— テンソルの核はイデアルの拡大 -/

/-- ★★★★**`ker (S → S ⊗_R (R⧸J)) = J·S`**。

★機構は `Algebra.TensorProduct.quotIdealMapEquivTensorQuot`
（`S ⧸ J·S ≃ₐ[S] S ⊗_R (R⧸J)`）1 本である。

★★これが「アフィンでの引き戻しが `Ideal.map` になる」ことの根拠である。 -/
theorem ker_algebraMap_tensor (R S : Type) [CommRing R] [CommRing S] [Algebra R S] (J : Ideal R) :
    RingHom.ker (algebraMap S (S ⊗[R] (R ⧸ J))) = J.map (algebraMap R S) := by
  set e := Algebra.TensorProduct.quotIdealMapEquivTensorQuot (A := R) S J with he
  have hcomp : (algebraMap S (S ⊗[R] (R ⧸ J)))
      = (e : (S ⧸ J.map (algebraMap R S)) →+* (S ⊗[R] (R ⧸ J))).comp
        (algebraMap S (S ⧸ J.map (algebraMap R S))) := by
    ext a; simp [he]
  have hker : RingHom.ker (e : (S ⧸ J.map (algebraMap R S)) →+* (S ⊗[R] (R ⧸ J))) = ⊥ := by
    ext x; simp
  rw [hcomp, ← RingHom.comap_ker, hker, Ideal.Quotient.algebraMap_eq,
    ← RingHom.ker_eq_comap_bot, Ideal.mk_ker]

/-! ## ★★★★★★★★アフィンでの `comap` -/

/-- ★★★★★★★★**閉部分スキームの引き戻しは `Spec(S ⊗_R R⧸J)` の核**。

★`ker_fst_of_isClosedImmersion`（mathlib）で `comap` を `pullback.fst` の核に直し、
`pullbackSpecIso_inv_fst′`（mathlib）でアフィンの引き戻しを
`Spec (S ⊗_R R⧸J)` に同定する。★★同型は `ker_comp_of_isIso` で落ちる。 -/
theorem comap_ker_quotient (R S : Type) [CommRing R] [CommRing S] [Algebra R S] (J : Ideal R) :
    (Scheme.Hom.ker (Spec.map (CommRingCat.ofHom (algebraMap R (R ⧸ J))))).comap
        (Spec.map (CommRingCat.ofHom (algebraMap R S)))
      = Scheme.Hom.ker (Spec.map (CommRingCat.ofHom (algebraMap S (S ⊗[R] (R ⧸ J))))) := by
  haveI : IsClosedImmersion (Spec.map (CommRingCat.ofHom (algebraMap R (R ⧸ J)))) := by
    rw [Ideal.Quotient.algebraMap_eq]
    exact IsClosedImmersion.spec_of_surjective _ Ideal.Quotient.mk_surjective
  have hfst : pullback.fst (Spec.map (CommRingCat.ofHom (algebraMap R S)))
        (Spec.map (CommRingCat.ofHom (algebraMap R (R ⧸ J))))
      = (pullbackSpecIso R S (R ⧸ J)).hom
        ≫ Spec.map (CommRingCat.ofHom (algebraMap S (S ⊗[R] (R ⧸ J)))) := by
    rw [← pullbackSpecIso_inv_fst' R S (R ⧸ J), Iso.hom_inv_id_assoc]
  rw [← Scheme.IdealSheafData.ker_fst_of_isClosedImmersion
    (Spec.map (CommRingCat.ofHom (algebraMap R (R ⧸ J))))
    (Spec.map (CommRingCat.ofHom (algebraMap R S))), hfst,
    Scheme.Hom.ker_comp_of_isIso]

/-- ★★★★★★★★★★**アフィンでの `ideal_comap`** ——
`mathlib-gap` の `ideal-comap-on-affine-opens` の**アフィンの場合**。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

> `(comap of V(J)).ideal ⊤ = J·S`（`ΓSpecIso` の分だけずれる）

★★**これが「引き戻した因子はイデアルの拡大である」の形式化**である。
★★★`comap_mul` はここから `Ideal.map_mul` で出る。 -/
theorem ideal_comap_ker_quotient_top (R S : Type) [CommRing R] [CommRing S] [Algebra R S]
    (J : Ideal R) :
    ((Scheme.Hom.ker (Spec.map (CommRingCat.ofHom (algebraMap R (R ⧸ J))))).comap
        (Spec.map (CommRingCat.ofHom (algebraMap R S)))).ideal
      ⟨⊤, AlgebraicGeometry.isAffineOpen_top (Spec (CommRingCat.of S))⟩
      = Ideal.comap (Scheme.ΓSpecIso (CommRingCat.of S)).hom.hom (J.map (algebraMap R S)) := by
  rw [comap_ker_quotient, ker_specMap_ideal_top]
  congr 1
  exact ker_algebraMap_tensor R S J

/-! ### ★出典の紐付け(`.src`)

★★これは `Remark 1.5.1` / `Remark 1.4.1` が要求する**配管**であって、
原典の命題そのものではない。 -/

def ker_specMap_ideal_top.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(配管——Spec.map の核をアフィンで読む)",
    sectionId := "genell-rem-1-5-1" }

def ker_algebraMap_tensor.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(配管——テンソルの核はイデアルの拡大)",
    sectionId := "genell-rem-1-5-1" }

def ideal_comap_ker_quotient_top.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(配管——アフィンでの ideal_comap。一般の場合は含まない)",
    sectionId := "genell-rem-1-5-1" }

def ideal_comap_ker_quotient_top.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Scheme.Hom.ker_apply(準コンパクトな射の核をアフィンで読む)"
      (.inMathlib "AlgebraicGeometry.Scheme.Hom.ker_apply") 9,
    .citation "[mathlib]"
      "Scheme.IdealSheafData.ker_fst_of_isClosedImmersion(comap は pullback.fst の核)"
      (.inMathlib "AlgebraicGeometry.Scheme.IdealSheafData.ker_fst_of_isClosedImmersion") 9,
    .citation "[mathlib]" "pullbackSpecIso_inv_fst′(アフィンの引き戻しはテンソル積)"
      (.inMathlib "AlgebraicGeometry.pullbackSpecIso_inv_fst'") 9,
    .citation "[mathlib]"
      "Algebra.TensorProduct.quotIdealMapEquivTensorQuot(S ⧸ J·S ≃ S ⊗_R R⧸J)"
      (.inMathlib "Algebra.TensorProduct.quotIdealMapEquivTensorQuot") 9,
    .implicitStep
      ("★★★★**部品はすべて mathlib にあった**——無かったのは繋ぎ方だけである。" ++
       "★mathlib には ideal_comap_of_isOpenImmersion(開埋め込みの場合)はあるが、" ++
       "アフィン射の場合も一般の場合も無い(2026-08-27 実測)") 9,
    .implicitStep
      ("★★残る段: 一般の(アフィンでない)X・Y の場合。" ++
       "★道筋は comap_comp で開埋め込みに沿って制限し、" ++
       "ideal_comap_of_isOpenImmersion で .ideal U に降ろし、" ++
       "本ファイルのアフィンの場合を当てる") 9 ]

end ABC3.Found.GenEll
