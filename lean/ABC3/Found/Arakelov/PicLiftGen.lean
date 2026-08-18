import ABC3.Found.Arakelov.PicLocGen

/-!
# Arakelov (B1) 第 130 ブロック —— **制限した生成元は生成元である**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★§9-145 の残る `sorry` の**中身そのもの**

§9-145 で「残るは **1 点**——制限した生成元が生成元であることの**同定**だけ」
と書いた、その 1 点である。

第 113 ブロック(`bijective_localizedGenerator`)は

    Bijective (IsLocalizedModule.map S' (Algebra.linearMap …) (liftAwayMapA …) e.symm)

を与えるが、★これは**抽象的な誘導射**であって「生成元による乗法」の形をしていない。

★★本ブロックはその **2 つが等しい**ことを言う:

    (c ↦ c • liftAwayMapA (e.symm 1))|_{R_g}  =  IsLocalizedModule.map … e.symm

★★★両辺は `R_g` 線型で `R_{t·g}` は `R_g` の局所化だから、
`IsLocalizedModule.ext` により **`algebraMap` に合成した後で一致**すればよい。

    左辺 (algebraMap c) = c • liftAwayMapA (e.symm 1) = liftAwayMapA (e.symm c)
    右辺 (algebraMap c) = liftAwayMapA (e.symm c)                      ✅

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `toSpanSingleton_eq_localizedMap` | ★★★乗法射 = 誘導射 |
| `bijective_smul_liftGen` | ★★★★★★**制限した生成元による乗法は全単射** |

## ★★★これで「基本開集合の族」で `IsLocallyTrivial` が言える

残るのは**切断側への輸送**(第 125・126 の同型で運ぶ)だけである。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite

variable (R : Type u) [CommRing R] (g t : R)
  (M : Type u) [AddCommGroup M] [Module R M]

/-- ★★★**「生成元による乗法」は「局所化した誘導射」に等しい**。 -/
theorem toSpanSingleton_eq_localizedMap
    (e : LocalizedModule (Submonoid.powers g) M
      ≃ₗ[Localization (Submonoid.powers g)] Localization (Submonoid.powers g)) :
    (LinearMap.toSpanSingleton (Localization (Submonoid.powers (t * g)))
        (LocalizedModule (Submonoid.powers (t * g)) M)
        (liftAwayMapA R g t M (e.symm 1))).restrictScalars
        (Localization (Submonoid.powers g))
      = IsLocalizedModule.map ((Submonoid.closure ({t, g} : Set R)).map
          (algebraMap R (Localization (Submonoid.powers g))))
        (Algebra.linearMap (Localization (Submonoid.powers g))
          (Localization (Submonoid.powers (t * g))))
        (liftAwayMapA R g t M) e.symm.toLinearMap := by
  refine IsLocalizedModule.ext ((Submonoid.closure ({t, g} : Set R)).map
      (algebraMap R (Localization (Submonoid.powers g))))
    (Algebra.linearMap (Localization (Submonoid.powers g))
      (Localization (Submonoid.powers (t * g))))
    (fun s => IsLocalizedModule.map_units (liftAwayMapA R g t M) s) ?_
  refine LinearMap.ext fun c => ?_
  show (algebraMap (Localization (Submonoid.powers g))
      (Localization (Submonoid.powers (t * g))) c) • liftAwayMapA R g t M (e.symm 1) = _
  rw [algebraMap_smul, ← map_smul, ← map_smul, smul_eq_mul, mul_one]
  exact (DFunLike.congr_fun (IsLocalizedModule.map_comp _ _
    (liftAwayMapA R g t M) e.symm.toLinearMap) c).symm

/-- ★★★★★★**制限した生成元による乗法は全単射である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが §9-145 で「残る 1 点」と書いたものである
——`D(g)` の上の生成元を `D(t·g)` へ制限してもやはり生成元である。 -/
theorem bijective_smul_liftGen
    (e : LocalizedModule (Submonoid.powers g) M
      ≃ₗ[Localization (Submonoid.powers g)] Localization (Submonoid.powers g)) :
    Function.Bijective (fun c : Localization (Submonoid.powers (t * g)) =>
      c • liftAwayMapA R g t M (e.symm 1)) := by
  have := bijective_localizedGenerator R g t M e
  rw [← toSpanSingleton_eq_localizedMap R g t M e] at this
  exact this

/-! ## ★出典の紐付け(`.src`) -/

def bijective_smul_liftGen.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——制限した生成元は生成元)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
