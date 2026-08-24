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

## ★本ファイルの内容(`Ω` の材料 —— 残るは `map_comp` 1 本)

| | 状態 |
|---|---|
| `biratRtIso` …… 根の食い違いを吸収する同型(＋三角形) | ★済 |
| `biratFT` / `biratDegEq` …… `homPfMap` の 2 入力 | ★済 |
| `omegaObj` / `omegaMap` …… `Ω` の対象と射 | ★済 |
| `omegaMap_id` …… `map_id` | ★済 |
| `biratRtIso_rtLift` …… `β` と根の持ち上げの四角形 | ★済 |
| `omegaMap_comp` …… `map_comp` | ☆残(道具は 4 つとも揃った) |

## ★★★★`omegaMap_comp` の手順書(2026-08-25 に測った)

★`Ω` の射を **`compPf` の共役ではなく `rootIso` の側**で書き直したので、
`map_comp` は**在庫 2 本＋新設 2 本**だけで組める。記号を
`Ψ = toBiratCat`、`hm = homPfMap F F' Ψ …`、`ρ(a,b) = rootIso a _ b _ _`、
`β_{A,d} = biratRtIso F F' A d` とする。

| 道具 | 出どころ |
|---|---|
| `homPfMap_rootIso_hom` …… `hm ∘ ρ.hom = ρ(Ψa,Ψb).hom ∘ hm` | ★本増分(`Thm34Pf.lean`) |
| `biratRtIso_rtLift` …… `rtLift^birat ≫ β_t = β_d ≫ Ψ(rtLift)` | ★本増分 |
| `rootIso_trans` …… `ρ(a).hom ∘ ρ(a').hom = ρ(a ≫ a').hom` | 在庫(`Def31Pf.lean`) |
| `homPfMap_compPf` / `rootIso_comp'` | 在庫 |

★★**鎖**(`compRoot f g = (rtRootIso …).hom (compPf (ρ⁻¹f) (ρ⁻¹g))` の外側から):

```
hm ((rtRootIso X Z).hom w)
  = ρ(Ψ(rtLift X), Ψ(rtLift Z)).hom (hm w)          -- homPfMap_rootIso_hom
ρ(β, β').hom (ρ(Ψ rtLift, Ψ rtLift').hom (hm w))
  = ρ(β ≫ Ψ(rtLift), β' ≫ Ψ(rtLift')).hom (hm w)    -- rootIso_trans
  = ρ(rtLift^birat ≫ β_t, rtLift'^birat ≫ β'_t).hom (hm w)   -- ★biratRtIso_rtLift
  = (rtRootIso^birat X Z).hom (ρ(β_t, β'_t).hom (hm w))       -- rootIso_trans(逆向き)
```

★内側は `homPfMap_compPf`(`hm` は `compPf` を保つ)と
`rootIso_comp'`(`ρ.hom` は `compPf` を分配する)で分けてから、
`.inv` の側は「`hm ∘ ρ.hom = ρ'.hom ∘ hm` ⟹ `hm ∘ ρ.inv = ρ'.inv ∘ hm`」
(両者が同型だから)で降ろす。

★★残る手間は**指数の帳簿**(`compRoot` の `mul_comm` 6 か所)だけである。
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
/-- ★★★**射** —— `homPfMap` で送って、根の食い違いを **`rootIso`** で吸収する。

★★★**設計の要点(2026-08-25 に取り直した)**: 最初は
`compPf` による共役(`β ≫ – ≫ β⁻¹`)で書いていたが、
**`rootIso`(根の不変性)の側で書くほうが良い** ——
`β` は同型なので Frobenius 型・次数 1 で、`rootIso` の 3 仮定をそのまま満たす。
★こう書くと `map_comp` が在庫の `rootIso_trans` / `rootIso_comp'` と
新設の `homPfMap_rootIso_hom` / `biratRtIso_rtLift` だけで組める。 -/
noncomputable def omegaMap (F' : FrobenioidCore (biratPre P G))
    {X Y : PfRootObj P F} (f : HomRoot P F X Y) :
    HomRoot (biratPre P G) F' (omegaObj F F' X) (omegaObj F F' Y) :=
  (rootIso (F := F') (biratRtIso F F' X.obj Y.root)
      (isFrobeniusType_of_isIso (biratPre P G) _)
      (biratRtIso F F' Y.obj X.root)
      (isFrobeniusType_of_isIso (biratPre P G) _)
      (by rw [isLinear_of_isIso (biratPre P G) _,
        isLinear_of_isIso (biratPre P G) _])).hom
    (homPfMap F F' (toBiratCat P G) biratFT biratDegEq _ _ f)

/-! ## ★4. `map_id` -/

variable (F) in
/-- ★★★**`Ω` は恒等射を保つ** —— `rootIso_hom_id` そのもの。 -/
theorem omegaMap_id (F' : FrobenioidCore (biratPre P G)) (X : PfRootObj P F) :
    omegaMap F F' (idRoot P F X) = idRoot (biratPre P G) F' (omegaObj F F' X) := by
  have hmap := homPfMap_toHomPf F F' (toBiratCat P G) biratFT biratDegEq
    (𝟙 (rtObj P F X.obj X.root))
  rw [(toBiratCat P G).map_id] at hmap
  refine Eq.trans (congrArg (rootIso (F := F') (biratRtIso F F' X.obj X.root)
      (isFrobeniusType_of_isIso (biratPre P G) _)
      (biratRtIso F F' X.obj X.root)
      (isFrobeniusType_of_isIso (biratPre P G) _)
      (by rw [isLinear_of_isIso (biratPre P G) _])).hom hmap) ?_
  exact rootIso_hom_id (F := F') (biratRtIso F F' X.obj X.root)
    (isFrobeniusType_of_isIso (biratPre P G) _)

/-! ## ★5. `map_comp` の要 —— **持ち上げとの両立** -/

variable (F) in
/-- ★★★★★★**`biratRtIso` は根の持ち上げと両立する**。

```
rtObj^birat(A,d) --rtLift^birat--> rtObj^birat(A,t)
      |β_d                              |β_t
      v                                 v
 (rtObj(A,d))^birat --Ψ(rtLift)--> (rtObj(A,t))^birat
```

★★**証明は 2 行の骨**である ——
`𝒞^birat` は totally epimorphic なので構造射 `rtExt` は epi、
両辺に前から `rtExt` を掛けると `rtLift_ext` と**三角形** `biratRtIso_tri` で
どちらも `Ψ(rtExt A t)` になる。

★★★これが `omegaMap_comp`(`rtLift` と `homPfMap` の可換性)の要である。
`β` が `Exists.choose` でも構わないのは、要るのが**この四角形だけ**だから。 -/
theorem biratRtIso_rtLift (F' : FrobenioidCore (biratPre P G)) (A : C)
    {d e t : ℕ+} (h : t = e * d) :
    rtLift (biratPre P G) F' (biratUp P G A) h ≫ biratRtIso F F' A t
      = biratRtIso F F' A d ≫ (toBiratCat P G).map (rtLift P F A h) := by
  haveI hepi : Epi (rtExt (biratPre P G) F' (biratUp P G A) d) :=
    (biratPre P G).totEpiC _ _ _
  refine (cancel_epi (rtExt (biratPre P G) F' (biratUp P G A) d)).mp ?_
  have hL : rtExt (biratPre P G) F' (biratUp P G A) d
        ≫ (rtLift (biratPre P G) F' (biratUp P G A) h ≫ biratRtIso F F' A t)
      = (toBiratCat P G).map (rtExt P F A t) := by
    rw [← Category.assoc, rtLift_ext, biratRtIso_tri]
  have hR : rtExt (biratPre P G) F' (biratUp P G A) d
        ≫ (biratRtIso F F' A d ≫ (toBiratCat P G).map (rtLift P F A h))
      = (toBiratCat P G).map (rtExt P F A t) := by
    have hc : (toBiratCat P G).map (rtExt P F A d) ≫ (toBiratCat P G).map (rtLift P F A h)
        = (toBiratCat P G).map (rtExt P F A d ≫ rtLift P F A h) :=
      ((toBiratCat P G).map_comp _ _).symm
    rw [← Category.assoc, biratRtIso_tri]
    exact hc.trans (congrArg (toBiratCat P G).map (rtLift_ext P F A h))
  exact hL.trans hR.symm

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
