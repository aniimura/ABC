import ABC3.Found.Arakelov.PicTensorFrac

/-!
# Arakelov (B1) 第 82 ブロック —— **切断のテンソル写像**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★前層射の `app`

第 81 ブロックで「点ごとのテンソルは局所分数条件を保つ」が出た。
★本ブロックはそれを **`𝒪(U)` 双線型写像**に組み、テンソル積からの線型写像にする:

    (tilde M)(U) ⊗_{𝒪(U)} (tilde N)(U)  ⟶  (tilde (M ⊗_R N))(U)

## ★★機構 —— 切断加群の構造は**点ごと**である(実測)

| 主張 | 結果 |
|---|---|
| `(s + s').1 x = s.1 x + s'.1 x` | ★**`rfl`** |
| `(c • s).1 x = c.1 x • s.1 x` | ★**`rfl`** |

★★したがって加法性・斉次性は**点ごとに落ちる**——
加法は `TensorProduct.add_tmul`、斉次性は `map_smul` + `smul_tmul'`。

## ★★本ブロックで取れるもの

| 定義 | 内容 |
|---|---|
| `tensorSection` | ★切断のテンソル |
| `tensorSectionBilin` | ★★`𝒪(U)` 双線型写像 |
| `tensorSectionMap` | ★★★★**テンソル積からの線型写像** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace TensorProduct StructureSheaf

variable (R : Type u) [CommRing R] (M N : Type u)
  [AddCommGroup M] [Module R M] [AddCommGroup N] [Module R N]
  (U : Opens (PrimeSpectrum.Top R))

/-! ## ★点ごとの性質 -/

theorem tensorFun_add_left (s s' : Π x : U, Localizations (R := R) M x.1)
    (t : Π x : U, Localizations (R := R) N x.1) :
    tensorFun R M N U (s + s') t = tensorFun R M N U s t + tensorFun R M N U s' t := by
  funext x; simp [tensorFun, TensorProduct.add_tmul]

theorem tensorFun_add_right (s : Π x : U, Localizations (R := R) M x.1)
    (t t' : Π x : U, Localizations (R := R) N x.1) :
    tensorFun R M N U s (t + t') = tensorFun R M N U s t + tensorFun R M N U s t' := by
  funext x; simp [tensorFun, TensorProduct.tmul_add]

theorem tensorFun_smul_left (c : Π x : U, Localizations (R := R) R x.1)
    (s : Π x : U, Localizations (R := R) M x.1)
    (t : Π x : U, Localizations (R := R) N x.1) (x : U) :
    tensorFun R M N U (fun y => c y • s y) t x = c x • tensorFun R M N U s t x := by
  simp only [tensorFun]; rw [← map_smul]; congr 1

theorem tensorFun_smul_right (c : Π x : U, Localizations (R := R) R x.1)
    (s : Π x : U, Localizations (R := R) M x.1)
    (t : Π x : U, Localizations (R := R) N x.1) (x : U) :
    tensorFun R M N U s (fun y => c y • t y) x = c x • tensorFun R M N U s t x := by
  simp only [tensorFun]; rw [← map_smul]; congr 1
  exact (TensorProduct.tmul_smul
    (R := Localization (x.1 : PrimeSpectrum R).asIdeal.primeCompl)
    (c x) (s x) (t x))

/-! ## ★★切断のテンソル -/

/-- ★**切断のテンソル**。 -/
noncomputable def tensorSection
    (s : (structureSheafInType R M).1.obj (op U))
    (t : (structureSheafInType R N).1.obj (op U)) :
    (structureSheafInType R (M ⊗[R] N)).1.obj (op U) :=
  ⟨tensorFun R M N U s.1 t.1, tensorFun_frac R M N U s.2 t.2⟩

/-- ★★★★**切断のテンソル写像**(`𝒪(U)` 双線型)。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが前層射の `app U` の中身である。 -/
noncomputable def tensorSectionBilin :
    (structureSheafInType R M).1.obj (op U)
      →ₗ[(structureSheafInType R R).1.obj (op U)]
    (structureSheafInType R N).1.obj (op U)
      →ₗ[(structureSheafInType R R).1.obj (op U)]
    (structureSheafInType R (M ⊗[R] N)).1.obj (op U) where
  toFun s :=
    { toFun := fun t => tensorSection R M N U s t
      map_add' := fun t t' => Subtype.ext (tensorFun_add_right R M N U s.1 t.1 t'.1)
      map_smul' := fun c t => Subtype.ext (funext fun x =>
        tensorFun_smul_right R M N U c.1 s.1 t.1 x) }
  map_add' s s' := LinearMap.ext fun t =>
    Subtype.ext (tensorFun_add_left R M N U s.1 s'.1 t.1)
  map_smul' c s := LinearMap.ext fun t =>
    Subtype.ext (funext fun x => tensorFun_smul_left R M N U c.1 s.1 t.1 x)

/-- ★★★★**テンソル積からの線型写像**。 -/
noncomputable def tensorSectionMap :
    TensorProduct ((structureSheafInType R R).1.obj (op U))
        ((structureSheafInType R M).1.obj (op U))
        ((structureSheafInType R N).1.obj (op U))
      →ₗ[(structureSheafInType R R).1.obj (op U)]
    (structureSheafInType R (M ⊗[R] N)).1.obj (op U) :=
  TensorProduct.lift (tensorSectionBilin R M N U)

/-! ## ★出典の紐付け(`.src`) -/

def tensorSectionMap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——切断のテンソル写像)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
