/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop55ScaleRootCoa
import ABC3.Found.FrdI.Thm34Pf

/-!
# [FrdI] Proposition 5.5, (ii) —— `Ω : 𝒞^pf ⥤ (𝒞^birat)^pf` の材料

原文 (FrdI p.104):
> (ii) There is a natural equivalence of categories [compatible with the functors to the

★★`Prop44Univ.lean` で **`𝒞^birat` の普遍性**(`biratDescFunctor`)を作った。
これに流すべき関手が

```
Ω : 𝒞^pf ⥤ (𝒞^birat)^pf
```

である。★`Ω` が co-angular pre-step を同型へ送ることさえ言えば、

```
Θ = biratDescFunctor Ω hΩ : (𝒞^pf)^birat ⥤ (𝒞^birat)^pf
```

が**関手として無料で出る**(合成の coherence を手で書かずに済む)。

## ★★測って分かったこと(2026-08-25)—— `frobDegUniq` は**三角形も返す**

`exists_rtObj_birat_iso`(`Prop55ScaleRootCoa.lean`)は
`F'.frobDegUniq` の返り値 `⟨β, hβ, -⟩` の**第 3 成分を捨てていた**。
★★捨てていたのは

```
rtExt (𝒞^birat) F' (A^birat) d ≫ β = (toBiratCat).map (rtExt 𝒞 F A d)
```

という**三角形**である。★`Ω` の `map_comp` は「根の持ち上げ `rtLift` と
`homPfMap` が可換」を要求し、その可換性はこの三角形から出る筋になっている。
★したがって `β` が `Exists.choose` であっても**構わない** ——
必要なのは `β` の一意性ではなく、この三角形だけである。

## ★本ファイルの内容(`Ω` の 4 材料のうち 3 つ)

| | 状態 |
|---|---|
| `biratFT` / `biratDegEq` …… `homPfMap` の 2 入力 | ★済 |
| `biratRtIso` …… 根の食い違いを吸収する同型(＋三角形) | ★済 |
| `omegaObj` / `omegaMap` …… `Ω` の対象と射 | ★済 |
| `omegaMap_id` …… `map_id` | ★済 |
| `omegaMap_comp` …… `map_comp`(`rtLift` と `homPfMap` の可換性) | ☆残 |
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section BiratOmega

/-! ### ★★★★測って分かったこと(2026-08-25)—— **universe の壁**と、その越え方

`homPfMap`(`Theorem 3.4, (iii)`、在庫)は
`C₁ C₂ : Type u2` を**同じ hom 宇宙 `v2`** で受ける。これは飾りではなく、
`HomPf` が `TypeCat.{max u2 v2}` への関手の**余極限**なので
**余錐の頂点が同じ宇宙になければならない**ためである。

★ところが `BiratCat P G` の hom は `Type (max v u2 v2)` にいる
(`biratCategory : Category.{max v u2 v2}`)。したがって
`max u2 v2 =?= max u2 v v2` が解けず、そのままでは `homPfMap` を流せない
(実測: `failed to solve universe constraint max u2 v2 =?= max v3 u3`)。

★★**越え方は 2 つある**:

1. `HomPf` の余極限からの写像を**より大きな宇宙へ**降ろす API を足す
   (`Types` の商としての記述＋`Quot.lift`)。★汎用だが `Def31Pf` の改修が要る。
2. ★**`𝒞` と `𝒟` の hom 宇宙を対象の宇宙 `u2` に揃える**(本ファイルはこちら)。
   ★`homPfMap` は `C₁` と `C₂` に**同じ `v2`** を要求し、
   `BiratCat` の hom は `max v u2 v2` なので、
   `v := v2 := u2` と取れば 3 つがすべて `u2` に潰れて条件が満たされる。
   これは数学的な仮定ではなく**宇宙パラメータの特殊化**であり、
   実際に使う Frobenioid は具体的な宇宙にいるので後続に影響しない。

★逸脱の記録: 本ファイルは `{C : Type u2} [Category.{u2} C]`,
`[Category.{u2} D]` を置く(原典の主張には何も足していない)。 -/

variable {D : Type u} [Category.{u2} D] {C : Type u2} [Category.{u2} C]
  {Φ : MonoidOn.{u2, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}
  {G : Frobenioid P} [IsConnected D]

/-! ## ★1. `homPfMap` の 2 入力 -/

/-- ★★**`toBiratCat` は Frobenius 型を Frobenius 型へ送る**。

★在庫 `birat_isFrobeniusType_iff` は「co-angular ∧ 底同型」に落とす。
`𝒞` の Frobenius 型はその 2 つを(等長性とともに)含んでいる。 -/
theorem biratFT {X Y : C} (f : X ⟶ Y) (h : IsFrobeniusType P f) :
    IsFrobeniusType (biratPre P G) ((toBiratCat P G).map f) :=
  (birat_isFrobeniusType_iff P G f).mpr ⟨h.1.1, h.2⟩

/-- ★**`toBiratCat` は Frobenius 次数を変えない**。 -/
theorem birat_degFr_map {X Y : C} (f : X ⟶ Y) :
    (biratPre P G).degFr ((toBiratCat P G).map f) = P.degFr f := by
  show biratDeg (toHomBirat (P := P) (G := G) f) = P.degFr f
  exact biratDeg_toHomBirat f

/-- ★★**`homPfMap` の第 2 入力** —— 次数の等号は保たれる。 -/
theorem biratDegEq {X Y X' Y' : C} (f : X ⟶ Y) (g : X' ⟶ Y')
    (h : P.degFr f = P.degFr g) :
    (biratPre P G).degFr ((toBiratCat P G).map f)
      = (biratPre P G).degFr ((toBiratCat P G).map g) := by
  rw [birat_degFr_map, birat_degFr_map, h]

/-! ## ★2. 根の食い違いを吸収する同型 —— **三角形つき** -/

/-- ★★★★★**`exists_rtObj_birat_iso` の強化版** ——
同型 `β` は**拡大射の三角形**も満たす。

★★`F'.frobDegUniq` の返り値は `∃ β, IsIso β ∧ φ ≫ β = ψ` で、
在庫はこの第 3 成分を捨てていた。★`Ω` の `map_comp` に要るのはこれである。 -/
theorem exists_rtObj_birat_iso_tri (F' : FrobenioidCore (biratPre P G))
    (Z : C) (d : ℕ+) :
    ∃ β : rtObj (biratPre P G) F' (biratUp P G Z) d ⟶ biratUp P G (rtObj P F Z d),
      IsIso β ∧ rtExt (biratPre P G) F' (biratUp P G Z) d ≫ β
        = (toBiratCat P G).map (rtExt P F Z d) := by
  have hfrob : IsFrobeniusType (biratPre P G) ((toBiratCat P G).map (rtExt P F Z d)) :=
    biratFT (G := G) _ (rtExt_frobType P F Z d)
  have hdeg : (biratPre P G).degFr (rtExt (biratPre P G) F' (biratUp P G Z) d)
      = (biratPre P G).degFr ((toBiratCat P G).map (rtExt P F Z d)) := by
    rw [rtExt_degFr, birat_degFr_map, rtExt_degFr]
  exact F'.frobDegUniq _ _ _
    (rtExt (biratPre P G) F' (biratUp P G Z) d)
    ((toBiratCat P G).map (rtExt P F Z d))
    (rtExt_frobType (biratPre P G) F' (biratUp P G Z) d) hfrob hdeg

variable (F) in
/-- ★★**選んだ同型** —— `(𝒞^birat)^pf` の根から `𝒞` の根の像へ。 -/
noncomputable def biratRtIso (F' : FrobenioidCore (biratPre P G)) (Z : C) (d : ℕ+) :
    rtObj (biratPre P G) F' (biratUp P G Z) d ⟶ biratUp P G (rtObj P F Z d) :=
  (exists_rtObj_birat_iso_tri (F := F) F' Z d).choose

variable (F) in
instance biratRtIso_isIso (F' : FrobenioidCore (biratPre P G)) (Z : C) (d : ℕ+) :
    IsIso (biratRtIso F F' Z d) :=
  (exists_rtObj_birat_iso_tri (F := F) F' Z d).choose_spec.1

variable (F) in
/-- ★★★**三角形** —— `Ω` の関手則で要る唯一の coherence。 -/
theorem biratRtIso_tri (F' : FrobenioidCore (biratPre P G)) (Z : C) (d : ℕ+) :
    rtExt (biratPre P G) F' (biratUp P G Z) d ≫ biratRtIso F F' Z d
      = (toBiratCat P G).map (rtExt P F Z d) :=
  (exists_rtObj_birat_iso_tri (F := F) F' Z d).choose_spec.2

/-! ## ★3. `Ω` の対象と射 -/

variable (F) in
/-- ★**対象** —— `⟨A, n⟩ ↦ ⟨A^birat, n⟩`。 -/
def omegaObj (F' : FrobenioidCore (biratPre P G)) (X : PfRootObj P F) :
    PfRootObj (biratPre P G) F' :=
  ⟨biratUp P G X.obj, X.root⟩

variable (F) in
/-- ★★★**射** —— `homPfMap` で送って、根の食い違いを `biratRtIso` で共役する。 -/
noncomputable def omegaMap (F' : FrobenioidCore (biratPre P G))
    {X Y : PfRootObj P F} (f : HomRoot P F X Y) :
    HomRoot (biratPre P G) F' (omegaObj F F' X) (omegaObj F F' Y) :=
  compPf (biratPre P G) F' (toHomPf (F := F') (biratRtIso F F' X.obj Y.root))
    (compPf (biratPre P G) F'
      (homPfMap F F' (toBiratCat P G) biratFT biratDegEq _ _ f)
      (toHomPf (F := F') (@inv _ _ _ _ (biratRtIso F F' Y.obj X.root) (biratRtIso_isIso F F' Y.obj X.root))))

/-! ## ★4. `map_id` -/

variable (F) in
/-- ★★★**`Ω` は恒等射を保つ** —— 共役が打ち消し合うだけ。

★★配管メモ: ここは `rw` で通らない。`rtObj` が `Exists.choose` なので
目標の中では `Classical.indefiniteDescription … |>.1` に展開されており、
`homPfMap … (rtObj P F X.obj X.root) …` という**畳んだ形のパターンが当たらない**。
★`congrArg` で項の側から組むと defeq で通る(idiom 14 と同じ逃げ方)。 -/
theorem omegaMap_id (F' : FrobenioidCore (biratPre P G)) (X : PfRootObj P F) :
    omegaMap F F' (idRoot P F X) = idRoot (biratPre P G) F' (omegaObj F F' X) := by
  have hmap := homPfMap_toHomPf F F' (toBiratCat P G) biratFT biratDegEq
    (𝟙 (rtObj P F X.obj X.root))
  rw [(toBiratCat P G).map_id] at hmap
  refine Eq.trans (congrArg (fun t => compPf (biratPre P G) F'
      (toHomPf (F := F') (biratRtIso F F' X.obj X.root))
      (compPf (biratPre P G) F' t
        (toHomPf (F := F') (@inv _ _ _ _ (biratRtIso F F' X.obj X.root)
          (biratRtIso_isIso F F' X.obj X.root))))) hmap) ?_
  refine Eq.trans (congrArg (compPf (biratPre P G) F'
      (toHomPf (F := F') (biratRtIso F F' X.obj X.root)))
    (compPf_id_left (P := biratPre P G) (F := F')
      (toHomPf (F := F') (@inv _ _ _ _ (biratRtIso F F' X.obj X.root)
        (biratRtIso_isIso F F' X.obj X.root))))) ?_
  refine Eq.trans (toHomPf_comp (F := F') (biratRtIso F F' X.obj X.root)
      (@inv _ _ _ _ (biratRtIso F F' X.obj X.root)
        (biratRtIso_isIso F F' X.obj X.root))).symm ?_
  exact congrArg (toHomPf (F := F')) (IsIso.hom_inv_id _)

end BiratOmega

/-! ### ★出典の紐付け -/

/-- ★★★★★locator —— `Proposition 5.5, (ii)` の `Ω : 𝒞^pf ⥤ (𝒞^birat)^pf` の
**根の食い違いを吸収する同型(三角形つき)**。 -/
def exists_rtObj_birat_iso_tri.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — 根の同型は拡大射の三角形も満たす",
    sectionId := "frdi-prop-5-5" }

/-- ★★★★★locator —— `Proposition 5.5, (ii)` の `Ω` の対象・射と `map_id`。 -/
def omegaMap_id.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — Ω : 𝒞^pf ⥤ (𝒞^birat)^pf の材料(map_id まで)",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
