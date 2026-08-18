import ABC3.Found.Arakelov.PicPrincipal
import ABC3.Found.Arakelov.PicEquivRing

/-!
# Arakelov (B2) 第 193 ブロック —— ★★★★★★**アフィンでは主因子＝単項イデアル**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★これで `CartierPicData` は 14 欄中 12 欄

    Cartier な D について:  idealSheaf D ≅ 𝒪_{Spec R}  ⟺  D.ideal ⊤ が単項

## ★★逆向きは B1 を丸ごと使った

★★★「大域切断が自由 ⟹ 層が自明」は**一般には偽**(準連接性が要る)。
本ブロックは B1 の `equivPicRingSheaf`(`Pic(Spec R) ≃* Pic R`、第 145)を使う:

    toPicRing R (ofDivisor D) = Pic.mk R (Γ(idealSheaf D)) = 1
      ⟹ ofDivisor D = 1  (同型なので単射)
      ⟹ idealSheaf D ≅ 𝒪  (第 192)

★B1 の中に既に「局所自明 ⟹ `fromTildeΓ` が同型」(第 143)が入っているので、
準連接性の議論を**もう一度書かずに済んだ**。

## ★★律速は係数環の綴りだった([[ring-instance-two-paths]])

`Γ(Spec R, ⊤)` と `R` は**同型だが等しくない**。`moduleSpecΓFunctor.obj` の
`R`-加群構造と、部分加群としての `Γ(⊤)`-加群構造を `restrictScalars` で繋ごうとすると
instance の経路が合わない。

★★**逃げ道**: 型付き恒等関数 `gammaVal` / `gammaValInv` を置くと

| 命題 | 結果 |
|---|---|
| `gammaVal (r • x) = (ΓSpecIso.inv r) • gammaVal x` | ★**`rfl`** |
| `gammaVal (x + y) = gammaVal x + gammaVal y` | ★**`rfl`** |
| 両側の合成が恒等 | ★**`rfl`** |

となり、線型同値を**手で組める**。第 173 と同じ手である。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `gammaVal` / `gammaValInv` | ★型の橋(両方向) |
| `gammaFreeEquiv` | ★★★大域切断加群は `R` 上自由 |
| `isPrincipalDivisor_affine` | ★★★★★★**アフィンでは主因子＝単項イデアル** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable (R : CommRingCat.{u}) (D : (Spec R).IdealSheafData)

def gammaVal (x : (AlgebraicGeometry.moduleSpecΓFunctor.obj (idealSheaf D) : Type u)) :
    ((idealSections D ⊤) : Type u) := x

theorem gammaVal_smul (r : (R : Type u))
    (x : (AlgebraicGeometry.moduleSpecΓFunctor.obj (idealSheaf D) : Type u)) :
    gammaVal R D (r • x) = ((Scheme.ΓSpecIso R).inv.hom r) • gammaVal R D x := rfl




def gammaValInv (x : ((idealSections D ⊤) : Type u)) :
    (AlgebraicGeometry.moduleSpecΓFunctor.obj (idealSheaf D) : Type u) := x

theorem gammaVal_add (x y : (AlgebraicGeometry.moduleSpecΓFunctor.obj (idealSheaf D) : Type u)) :
    gammaVal R D (x + y) = gammaVal R D x + gammaVal R D y := rfl

theorem gammaVal_gammaValInv (x : ((idealSections D ⊤) : Type u)) :
    gammaVal R D (gammaValInv R D x) = x := rfl

theorem gammaValInv_gammaVal (x : (AlgebraicGeometry.moduleSpecΓFunctor.obj (idealSheaf D) : Type u)) :
    gammaValInv R D (gammaVal R D x) = x := rfl




/-- ★★大域切断加群は `R` 上自由(生成元があるとき)。 -/
noncomputable def gammaFreeEquiv
    (eG : (Γ(Spec R, ⊤) : Type u) ≃ₗ[(Γ(Spec R, ⊤) : Type u)] ((idealSections D ⊤) : Type u)) :
    (AlgebraicGeometry.moduleSpecΓFunctor.obj (idealSheaf D) : Type u)
      ≃ₗ[(R : Type u)] (R : Type u) where
  toFun x := (Scheme.ΓSpecIso R).hom.hom (eG.symm (gammaVal R D x))
  map_add' a b := by
    rw [gammaVal_add, map_add, map_add]
  map_smul' r x := by
    show (Scheme.ΓSpecIso R).hom.hom (eG.symm (gammaVal R D (r • x))) = r * _
    rw [gammaVal_smul, map_smul]
    show (Scheme.ΓSpecIso R).hom.hom (((Scheme.ΓSpecIso R).inv.hom r) * _) = r * _
    rw [map_mul]
    congr 1
    exact congrArg (fun (m : _ ⟶ _) => (CommRingCat.Hom.hom m) r) (Scheme.ΓSpecIso R).inv_hom_id
  invFun r := gammaValInv R D (eG ((Scheme.ΓSpecIso R).inv.hom r))
  left_inv x := by
    have h1 : (Scheme.ΓSpecIso R).inv.hom ((Scheme.ΓSpecIso R).hom.hom (eG.symm (gammaVal R D x)))
        = eG.symm (gammaVal R D x) :=
      congrArg (fun (m : _ ⟶ _) => (CommRingCat.Hom.hom m) (eG.symm (gammaVal R D x)))
        (Scheme.ΓSpecIso R).hom_inv_id
    show gammaValInv R D (eG ((Scheme.ΓSpecIso R).inv.hom
      ((Scheme.ΓSpecIso R).hom.hom (eG.symm (gammaVal R D x))))) = x
    rw [h1, eG.apply_symm_apply]
    exact gammaValInv_gammaVal R D x
  right_inv r := by
    show (Scheme.ΓSpecIso R).hom.hom (eG.symm (gammaVal R D
      (gammaValInv R D (eG ((Scheme.ΓSpecIso R).inv.hom r))))) = r
    rw [gammaVal_gammaValInv, eG.symm_apply_apply]
    exact congrArg (fun (m : _ ⟶ _) => (CommRingCat.Hom.hom m) r) (Scheme.ΓSpecIso R).inv_hom_id



theorem isPrincipalDivisor_affine (R : CommRingCat.{u}) (D : (Spec R).IdealSheafData)
    (hD : IsCartier (Spec R) D) :
    IsPrincipalDivisor (Spec R) D ↔
      (D.ideal ⟨⊤, isAffineOpen_top (Spec R)⟩).IsPrincipal := by
  classical
  constructor
  · rintro ⟨e⟩
    have hli := (modulesIsoApp e (op (⊤ : (Spec R).Opens))).toLinearEquiv
    have hp := principal_of_equiv (idealSections D ⊤) hli.symm
    rw [idealSections_affine D ⟨⊤, isAffineOpen_top (Spec R)⟩] at hp
    exact hp
  · intro hp
    haveI := invertible_idealSections D ⟨⊤, isAffineOpen_top (Spec R)⟩ hD
    have hp' : (idealSections D ⊤).IsPrincipal := by
      rw [idealSections_affine D ⟨⊤, isAffineOpen_top (Spec R)⟩]
      exact hp
    obtain ⟨eG⟩ := invertible_principal_equiv (idealSections D ⊤) hp'
    haveI hinv : Module.Invertible (R : Type u)
        (AlgebraicGeometry.moduleSpecΓFunctor.obj (idealSheaf D) : Type u) :=
      invertible_gammaCarrier R (InvSheaf.ofLocallyTrivial (idealSheaf D)
        (isLocallyTrivial_idealSheaf D hD))
    have hod : ofDivisorSheaf D
        = PicSheaf.mk (InvSheaf.ofLocallyTrivial (idealSheaf D)
            (isLocallyTrivial_idealSheaf D hD)) := by
      rw [ofDivisorSheaf, divisorInvSheaf, dif_pos hD]
    have h1 : toPicRing R (ofDivisorSheaf D) = 1 := by
      rw [hod]
      show CommRing.Pic.mk (R : Type u)
        (AlgebraicGeometry.moduleSpecΓFunctor.obj (idealSheaf D) : Type u) = 1
      exact CommRing.Pic.mk_eq_one_iff.2 ⟨gammaFreeEquiv R D eG⟩
    have h2 : ofDivisorSheaf D = 1 := by
      refine (equivPicRingSheaf R).injective ?_
      rw [map_one]
      exact h1
    exact (ofDivisorSheaf_eq_one_iff D hD).1 h2


/-! ## ★出典の紐付け(`.src`) -/

def isPrincipalDivisor_affine.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——アフィンでは主因子は単項イデアル)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
