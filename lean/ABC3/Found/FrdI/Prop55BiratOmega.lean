/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop55ScaleRootCoa
import ABC3.Found.FrdI.Thm34Pf
import ABC3.Found.FrdI.Prop55PfArbFull
import ABC3.Found.FrdI.Prop44Univ

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

## ★本ファイルの内容 —— **`Ω` は関手として完成した**

| | 状態 |
|---|---|
| `biratRtIso` …… 根の食い違いを吸収する同型(＋三角形) | ★済 |
| `biratFT` / `biratDegEq` …… `homPfMap` の 2 入力 | ★済 |
| `omegaObj` / `omegaMap` …… `Ω` の対象と射 | ★済 |
| `omegaMap_id` …… `map_id` | ★済 |
| `biratRtIso_rtLift` …… `β` と根の持ち上げの四角形 | ★済 |
| `omegaMap_comp` / `omegaFunctor` …… `map_comp` と関手 | ★済 |

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

## ★★★★★残るのは `Θ.map` の全単射性 —— 測ったこと(2026-08-25)

`Θ` は組み上がり、**対象の上の全射性は無料**(`thetaFunctor_obj_surjective`)。
残るのは `Θ.map` が全単射(＝充満忠実)であること。

### ★★測定 1 —— **`Ω` から継承する道は無い**

* `Ω` は**充満でない**。`Hom_{(𝒞^pf)^birat}` は `Hom_{𝒞^pf}` より真に大きく、
  `(𝒞^birat)^pf` の射はその大きいほうと対応する(それが本命題の主張である)。
* `Ω` は**忠実でもない**。`toBiratCat` の射の写像は局所化だから単射でない。

★したがって `Θ` の充満忠実性は **在庫 `biratPfHom_bijective`** から来るしかない。

### ★★★測定 2 —— **要る等式は 1 本だけ**

`equiv := biratPfHomEquivRootFull : Hom_{(𝒞^birat)^pf}(ΩX,ΩY) ≃ Hom_{(𝒞^pf)^birat}(X,Y)`
は**全単射**(在庫、仮定なし、一般の根)。`Θ.map` はその逆向きの写像である。

★★**`∀ y, Θ.map (equiv y) = y` の 1 本だけ**を示せばよい ——
`equiv` が全射なので、これだけで `Θ.map = equiv.symm` が出て**全単射**になる。
★単射性・全射性を別々に論じる必要は**ない**。

### ★橋渡しの形(根 1 の段)

`biratPfMk hfi Gpf W Z ψ = HomBirat.mk (biratPfIdx hfi Gpf W Z)
  (rootMap hfi ψ (biratPfDeg W) ≫ (biratPfIsoB hfi W).inv)` なので、要るのは

```
Θ.map (biratPfMk hfi Gpf W Z ψ) = HomPf.mk ((idxToBirat …).obj W) (HomBirat.mk Z ψ)
```

★左辺は `biratDescHom_mk` で `inv (Ω (biratPfIdx …).hom) ≫ Ω (rootMap … ≫ isoB⁻¹)`。
★したがって新たに要るのは **`Ω` が `rootMap` / `rootLift` をどう写すか**の計算則で、
`omegaMap_mk`(代表元での表示)がその入口である。
残る 3 段(`scaleRootBiratHomEquiv` など)は `Θ` が関手であることから従う。
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
/-- ★★★★**`β` が定める根の不変性の同型**。

`β` は同型なので Frobenius 型・次数 1 で、`rootIso` の 3 仮定をそのまま満たす。 -/
noncomputable def betaIso (F' : FrobenioidCore (biratPre P G)) (A B : C) (dA dB : ℕ+) :
    HomPf (biratPre P G) F' (biratUp P G (rtObj P F A dA))
        (biratUp P G (rtObj P F B dB))
      ≅ HomPf (biratPre P G) F' (rtObj (biratPre P G) F' (biratUp P G A) dA)
        (rtObj (biratPre P G) F' (biratUp P G B) dB) :=
  rootIso (F := F') (biratRtIso F F' A dA) (isFrobeniusType_of_isIso (biratPre P G) _)
    (biratRtIso F F' B dB) (isFrobeniusType_of_isIso (biratPre P G) _)
    (by rw [isLinear_of_isIso (biratPre P G) _, isLinear_of_isIso (biratPre P G) _])

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
  (betaIso F F' X.obj Y.obj Y.root X.root).hom
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

/-! ## ★6. `Ω` と根の不変性の交換則 —— `map_comp` の全体 -/

/-- ★射が等しければ `rootIso` の `hom` も等しい(仮定は `Prop` なので `subst` で消える)。 -/
theorem rootIso_hom_congr (F' : FrobenioidCore (biratPre P G))
    {A B A' B' : BiratCat P G} {a a' : A ⟶ A'} (haa : a = a')
    (ha : IsFrobeniusType (biratPre P G) a) (ha' : IsFrobeniusType (biratPre P G) a')
    {b b' : B ⟶ B'} (hbb : b = b')
    (hb : IsFrobeniusType (biratPre P G) b) (hb' : IsFrobeniusType (biratPre P G) b')
    (hd : (biratPre P G).degFr a = (biratPre P G).degFr b)
    (hd' : (biratPre P G).degFr a' = (biratPre P G).degFr b')
    (z : HomPf (biratPre P G) F' A' B') :
    (rootIso (F := F') a ha b hb hd).hom z = (rootIso (F := F') a' ha' b' hb' hd').hom z := by
  subst haa
  subst hbb
  rfl

variable (F) in
/-- ★★★★★★★**`Ω` は根の不変性と交換する** ——

```
βIso(dA,dB).hom ∘ hm ∘ (rtRootIso 𝒞).hom
  = (rtRootIso 𝒞^birat).hom ∘ βIso(tA,tB).hom ∘ hm
```

★★**筋は 4 段**(ファイル冒頭の手順書のとおり):
1. `homPfMap_rootIso_hom` で `hm` を `rootIso` の外へ出す
2. `rootIso_trans` で 2 段を 1 段にする
3. ★**四角形 `biratRtIso_rtLift`** で `(β ≫ Ψ rtLift)` を `(rtLift^birat ≫ β)` に替える
4. `rootIso_trans` を逆向きに使って 1 段を 2 段に戻す -/
theorem betaIso_rtRootIso (F' : FrobenioidCore (biratPre P G)) (A B : C)
    {dA dB e tA tB : ℕ+} (hA : tA = e * dA) (hB : tB = e * dB)
    (x : HomPf P F (rtObj P F A tA) (rtObj P F B tB)) :
    (betaIso F F' A B dA dB).hom
        (homPfMap F F' (toBiratCat P G) biratFT biratDegEq _ _
          ((rtRootIso P F A B hA hB).hom x))
      = (rtRootIso (biratPre P G) F' (biratUp P G A) (biratUp P G B) hA hB).hom
          ((betaIso F F' A B tA tB).hom
            (homPfMap F F' (toBiratCat P G) biratFT biratDegEq _ _ x)) := by
  -- ★4 種の Frobenius 性
  have hfA : IsFrobeniusType (biratPre P G) ((toBiratCat P G).map (rtLift P F A hA)) :=
    biratFT (G := G) _ (rtLift_frobType P F A hA)
  have hfB : IsFrobeniusType (biratPre P G) ((toBiratCat P G).map (rtLift P F B hB)) :=
    biratFT (G := G) _ (rtLift_frobType P F B hB)
  have hbA : IsFrobeniusType (biratPre P G) (biratRtIso F F' A dA) :=
    isFrobeniusType_of_isIso (biratPre P G) _
  have hbB : IsFrobeniusType (biratPre P G) (biratRtIso F F' B dB) :=
    isFrobeniusType_of_isIso (biratPre P G) _
  have hbA' : IsFrobeniusType (biratPre P G) (biratRtIso F F' A tA) :=
    isFrobeniusType_of_isIso (biratPre P G) _
  have hbB' : IsFrobeniusType (biratPre P G) (biratRtIso F F' B tB) :=
    isFrobeniusType_of_isIso (biratPre P G) _
  have hLA : IsFrobeniusType (biratPre P G)
      (rtLift (biratPre P G) F' (biratUp P G A) hA) := rtLift_frobType _ _ _ hA
  have hLB : IsFrobeniusType (biratPre P G)
      (rtLift (biratPre P G) F' (biratUp P G B) hB) := rtLift_frobType _ _ _ hB
  -- ★次数の一致
  have hdb : (biratPre P G).degFr (biratRtIso F F' A dA)
      = (biratPre P G).degFr (biratRtIso F F' B dB) := by
    rw [isLinear_of_isIso (biratPre P G) _, isLinear_of_isIso (biratPre P G) _]
  have hdb' : (biratPre P G).degFr (biratRtIso F F' A tA)
      = (biratPre P G).degFr (biratRtIso F F' B tB) := by
    rw [isLinear_of_isIso (biratPre P G) _, isLinear_of_isIso (biratPre P G) _]
  have hdl : P.degFr (rtLift P F A hA) = P.degFr (rtLift P F B hB) := by
    rw [rtLift_degFr, rtLift_degFr]
  have hdf : (biratPre P G).degFr ((toBiratCat P G).map (rtLift P F A hA))
      = (biratPre P G).degFr ((toBiratCat P G).map (rtLift P F B hB)) :=
    biratDegEq (G := G) (rtLift P F A hA) (rtLift P F B hB) hdl
  have hdL : (biratPre P G).degFr (rtLift (biratPre P G) F' (biratUp P G A) hA)
      = (biratPre P G).degFr (rtLift (biratPre P G) F' (biratUp P G B) hB) := by
    rw [rtLift_degFr, rtLift_degFr]
  -- ★合成の Frobenius 性と次数
  have hcA : IsFrobeniusType (biratPre P G)
      (biratRtIso F F' A dA ≫ (toBiratCat P G).map (rtLift P F A hA)) :=
    IsFrobeniusType.comp (biratPre P G) F' hbA hfA
  have hcB : IsFrobeniusType (biratPre P G)
      (biratRtIso F F' B dB ≫ (toBiratCat P G).map (rtLift P F B hB)) :=
    IsFrobeniusType.comp (biratPre P G) F' hbB hfB
  have hcA' : IsFrobeniusType (biratPre P G)
      (rtLift (biratPre P G) F' (biratUp P G A) hA ≫ biratRtIso F F' A tA) :=
    IsFrobeniusType.comp (biratPre P G) F' hLA hbA'
  have hcB' : IsFrobeniusType (biratPre P G)
      (rtLift (biratPre P G) F' (biratUp P G B) hB ≫ biratRtIso F F' B tB) :=
    IsFrobeniusType.comp (biratPre P G) F' hLB hbB'
  have hdc : (biratPre P G).degFr
        (biratRtIso F F' A dA ≫ (toBiratCat P G).map (rtLift P F A hA))
      = (biratPre P G).degFr
        (biratRtIso F F' B dB ≫ (toBiratCat P G).map (rtLift P F B hB)) := by
    rw [(biratPre P G).degFr_comp, (biratPre P G).degFr_comp, hdb]
    exact congrArg (fun t => t * (biratPre P G).degFr (biratRtIso F F' B dB)) hdf
  have hdc' : (biratPre P G).degFr
        (rtLift (biratPre P G) F' (biratUp P G A) hA ≫ biratRtIso F F' A tA)
      = (biratPre P G).degFr
        (rtLift (biratPre P G) F' (biratUp P G B) hB ≫ biratRtIso F F' B tB) := by
    rw [(biratPre P G).degFr_comp, (biratPre P G).degFr_comp, hdb', hdL]
  -- ★段 1 —— `hm` を `rootIso` の外へ出す
  have hstep1 := homPfMap_rootIso_hom F F' (toBiratCat P G) biratFT biratDegEq
    (rtLift P F A hA) (rtLift_frobType P F A hA)
    (rtLift P F B hB) (rtLift_frobType P F B hB) hdl x
  -- ★段 2 —— 2 段を 1 段に
  have step12 : (betaIso F F' A B dA dB).hom
        (homPfMap F F' (toBiratCat P G) biratFT biratDegEq _ _
          ((rtRootIso P F A B hA hB).hom x))
      = (rootIso (F := F')
          (biratRtIso F F' A dA ≫ (toBiratCat P G).map (rtLift P F A hA)) hcA
          (biratRtIso F F' B dB ≫ (toBiratCat P G).map (rtLift P F B hB)) hcB hdc).hom
        (homPfMap F F' (toBiratCat P G) biratFT biratDegEq _ _ x) := by
    refine Eq.trans (congrArg (fun t => (betaIso F F' A B dA dB).hom t) hstep1) ?_
    exact rootIso_trans (F := F') (biratRtIso F F' A dA) hbA (biratRtIso F F' B dB) hbB hdb
      ((toBiratCat P G).map (rtLift P F A hA)) hfA
      ((toBiratCat P G).map (rtLift P F B hB)) hfB hdf hcA hcB hdc _
  -- ★段 4 —— 右辺も 1 段に
  have step4 : (rtRootIso (biratPre P G) F' (biratUp P G A) (biratUp P G B) hA hB).hom
        ((betaIso F F' A B tA tB).hom
          (homPfMap F F' (toBiratCat P G) biratFT biratDegEq _ _ x))
      = (rootIso (F := F')
          (rtLift (biratPre P G) F' (biratUp P G A) hA ≫ biratRtIso F F' A tA) hcA'
          (rtLift (biratPre P G) F' (biratUp P G B) hB ≫ biratRtIso F F' B tB) hcB' hdc').hom
        (homPfMap F F' (toBiratCat P G) biratFT biratDegEq _ _ x) :=
    rootIso_trans (F := F') (rtLift (biratPre P G) F' (biratUp P G A) hA) hLA
      (rtLift (biratPre P G) F' (biratUp P G B) hB) hLB hdL
      (biratRtIso F F' A tA) hbA' (biratRtIso F F' B tB) hbB' hdb' hcA' hcB' hdc' _
  -- ★段 3 —— 四角形で 2 つの 1 段を同一視
  rw [step12, step4]
  exact rootIso_hom_congr F' (biratRtIso_rtLift F F' A hA).symm hcA hcA'
    (biratRtIso_rtLift F F' B hB).symm hcB hcB' hdc hdc' _

variable (F) in
/-- ★★同上、`inv` の側(同型なので両辺に `hom` を当てるだけ)。 -/
theorem betaIso_rtRootIso_inv (F' : FrobenioidCore (biratPre P G)) (A B : C)
    {dA dB e tA tB : ℕ+} (hA : tA = e * dA) (hB : tB = e * dB)
    (y : HomPf P F (rtObj P F A dA) (rtObj P F B dB)) :
    (betaIso F F' A B tA tB).hom
        (homPfMap F F' (toBiratCat P G) biratFT biratDegEq _ _
          ((rtRootIso P F A B hA hB).inv y))
      = (rtRootIso (biratPre P G) F' (biratUp P G A) (biratUp P G B) hA hB).inv
          ((betaIso F F' A B dA dB).hom
            (homPfMap F F' (toBiratCat P G) biratFT biratDegEq _ _ y)) := by
  have h := betaIso_rtRootIso F F' A B hA hB ((rtRootIso P F A B hA hB).inv y)
  rw [Iso.inv_hom_id_apply] at h
  rw [h, Iso.hom_inv_id_apply]

/-! ## ★7. `map_comp` -/

variable (F) in
/-- ★★★★★★★**`Ω` は合成を保つ**(`map_comp`)。

★★**筋**(ファイル冒頭の手順書のとおり):
* **外側** …… `betaIso_rtRootIso` で `(rtRootIso 𝒞).hom` を `(rtRootIso 𝒞^birat).hom` に替える
* **内側** …… `homPfMap_compPf`(`hm` は `compPf` を保つ)と
  `rootIso_comp'`(`ρ.hom` は `compPf` に分配する)で 2 本に分け、
  各本を `betaIso_rtRootIso_inv` で `(rtRootIso 𝒞^birat).inv (Ω –)` に替える

★指数の帳簿は `compRoot` のものをそのまま写すだけ(`Prop` なので `rfl` / `mul_comm`)。 -/
theorem omegaMap_comp (F' : FrobenioidCore (biratPre P G)) {X Y Z : PfRootObj P F}
    (f : HomRoot P F X Y) (g : HomRoot P F Y Z) :
    omegaMap F F' (compRoot P F f g)
      = compRoot (biratPre P G) F' (omegaMap F F' f) (omegaMap F F' g) := by
  have hu := betaIso_rtRootIso_inv F F' X.obj Y.obj
    (show Z.root * Y.root = Z.root * Y.root from rfl)
    (show Z.root * X.root = Z.root * X.root from rfl) f
  have hv := betaIso_rtRootIso_inv F F' Y.obj Z.obj
    (show Z.root * X.root = X.root * Z.root from mul_comm _ _)
    (show Y.root * X.root = X.root * Y.root from mul_comm _ _) g
  have hout := betaIso_rtRootIso F F' X.obj Z.obj
    (show Z.root * Y.root = Y.root * Z.root from mul_comm _ _)
    (show Y.root * X.root = Y.root * X.root from rfl)
    (compPf P F
      ((rtRootIso P F X.obj Y.obj
        (show Z.root * Y.root = Z.root * Y.root from rfl)
        (show Z.root * X.root = Z.root * X.root from rfl)).inv f)
      ((rtRootIso P F Y.obj Z.obj
        (show Z.root * X.root = X.root * Z.root from mul_comm _ _)
        (show Y.root * X.root = X.root * Y.root from mul_comm _ _)).inv g))
  refine Eq.trans hout ?_
  refine congrArg (fun t => (rtRootIso (biratPre P G) F' (biratUp P G X.obj)
      (biratUp P G Z.obj)
      (show Z.root * Y.root = Y.root * Z.root from mul_comm _ _)
      (show Y.root * X.root = Y.root * X.root from rfl)).hom t) ?_
  refine Eq.trans (congrArg
    (fun t => (betaIso F F' X.obj Z.obj (Z.root * Y.root) (Y.root * X.root)).hom t)
    (homPfMap_compPf F F' (toBiratCat P G) biratFT biratDegEq _ _)) ?_
  refine Eq.trans (rootIso_comp' (F := F')
      (biratRtIso F F' X.obj (Z.root * Y.root))
      (isFrobeniusType_of_isIso (biratPre P G) _)
      (biratRtIso F F' Y.obj (Z.root * X.root))
      (isFrobeniusType_of_isIso (biratPre P G) _)
      (biratRtIso F F' Z.obj (Y.root * X.root))
      (isFrobeniusType_of_isIso (biratPre P G) _)
      (by rw [isLinear_of_isIso (biratPre P G) _, isLinear_of_isIso (biratPre P G) _])
      (by rw [isLinear_of_isIso (biratPre P G) _, isLinear_of_isIso (biratPre P G) _])
      _ _).symm ?_
  exact congr (congrArg (compPf (biratPre P G) F') hu) hv


/-! ## ★8. 関手 `Ω` -/

variable (F) in
/-- ★★★★★★★**[FrdI] Proposition 5.5, (ii)** —— **`Ω : 𝒞^pf ⥤ (𝒞^birat)^pf`**。

★対象は `⟨A,n⟩ ↦ ⟨A^birat, n⟩`、射は `homPfMap` を `betaIso` で共役したもの。
★★これを `Prop44Univ.lean` の `biratDescFunctor`(`𝒞^birat` の普遍性)に流せば
`Θ : (𝒞^pf)^birat ⥤ (𝒞^birat)^pf` が**関手として出る**。 -/
noncomputable def omegaFunctor (F' : FrobenioidCore (biratPre P G)) :
    PfRootObj P F ⥤ PfRootObj (biratPre P G) F' where
  obj X := omegaObj F F' X
  map f := omegaMap F F' f
  map_id X := omegaMap_id F F' X
  map_comp f g := omegaMap_comp F F' f g

@[simp] theorem omegaFunctor_obj (F' : FrobenioidCore (biratPre P G)) (X : PfRootObj P F) :
    (omegaFunctor F F').obj X = omegaObj F F' X := rfl

@[simp] theorem omegaFunctor_map (F' : FrobenioidCore (biratPre P G))
    {X Y : PfRootObj P F} (f : HomRoot P F X Y) :
    (omegaFunctor F F').map f = omegaMap F F' f := rfl


/-! ## ★9. `Ω` は co-angular pre-step を同型に送る

★★★**測って分かったこと(2026-08-25)**: **co-angularity は使わない**。
`(𝒞^birat)^pf` は

* **すべての射が isometric**(`𝒞^birat` の零因子の単系が `PUnit` だから)
* **isotropic 型**(`𝒞` が isotropic なら `𝒞^birat` もそう、
  さらに `pfRoot_isOfIsotropicType` で `(𝒞^birat)^pf` もそう)

なので、**pre-step であるだけで同型になる**。
★したがって要るのは「`Ω` は pre-step を pre-step に送る」だけで、
それは**代表元の判定**(`isPreStep_mk_iff`)で `𝒞` の側に落ちる。 -/

variable (F) in
/-- ★★★**`Ω` の代表元での表示** —— 射は `Ψ` を当てるだけ、添字は 2 段押し出す。 -/
theorem omegaMap_mk (F' : FrobenioidCore (biratPre P G)) {X Y : PfRootObj P F}
    (Z : IdxPf P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root))
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) :
    omegaMap F F' (show HomRoot P F X Y from HomPf.mk Z φ)
      = HomPf.mk ((pushIdx (F := F') (biratRtIso F F' X.obj Y.root)
            (isFrobeniusType_of_isIso (biratPre P G) _)
            (biratRtIso F F' Y.obj X.root)
            (isFrobeniusType_of_isIso (biratPre P G) _)
            (by rw [isLinear_of_isIso (biratPre P G) _,
              isLinear_of_isIso (biratPre P G) _])).obj
          ((idxPfMap F F' (toBiratCat P G) biratFT biratDegEq _ _).obj Z))
        ((toBiratCat P G).map φ) := by
  show (betaIso F F' X.obj Y.obj Y.root X.root).hom
      (homPfMap F F' (toBiratCat P G) biratFT biratDegEq _ _ (HomPf.mk Z φ)) = _
  rw [homPfMap_mk]
  exact rootIso_hom_mk _ _ _ _ _ _ _

variable (F) in
/-- ★★★★★★**`Ω` は根 1 の添字の射を `β ≫ Ψ(–) ≫ β⁻¹` に送る**。

★★これが `pfKappa` / `toRootHom` / `rootMap` の像を計算する入口である
(`Θ.map` の全単射性の橋渡しで要る)。

★★配管メモ: 添字**対象の等式**では `rw` できないので、
`homPfMap_toHomPf` と同じく**添字圏の射 `(β,β′)` を 1 本作って
`HomPf.mk_map` で移す**。 -/
theorem omegaMap_toHomPf (F' : FrobenioidCore (biratPre P G)) {X Y : PfRootObj P F}
    (χ : rtObj P F X.obj Y.root ⟶ rtObj P F Y.obj X.root) :
    omegaMap F F' (show HomRoot P F X Y from toHomPf (F := F) χ)
      = (show HomRoot (biratPre P G) F' (omegaObj F F' X) (omegaObj F F' Y) from
          toHomPf (F := F') (biratRtIso F F' X.obj Y.root ≫ (toBiratCat P G).map χ
            ≫ @inv _ _ _ _ (biratRtIso F F' Y.obj X.root)
              (biratRtIso_isIso F F' Y.obj X.root))) := by
  have hbA : IsFrobeniusType (biratPre P G) (biratRtIso F F' X.obj Y.root) :=
    isFrobeniusType_of_isIso (biratPre P G) _
  have hbB : IsFrobeniusType (biratPre P G) (biratRtIso F F' Y.obj X.root) :=
    isFrobeniusType_of_isIso (biratPre P G) _
  have hdb : (biratPre P G).degFr (biratRtIso F F' X.obj Y.root)
      = (biratPre P G).degFr (biratRtIso F F' Y.obj X.root) := by
    rw [isLinear_of_isIso (biratPre P G) _, isLinear_of_isIso (biratPre P G) _]
  -- ★添字圏の射 `(β, β′)`
  let hv : idxOne (biratPre P G) F' (rtObj (biratPre P G) F' (biratUp P G X.obj) Y.root)
        (rtObj (biratPre P G) F' (biratUp P G Y.obj) X.root)
      ⟶ (pushIdx (F := F') (biratRtIso F F' X.obj Y.root) hbA
          (biratRtIso F F' Y.obj X.root) hbB hdb).obj
        ((idxPfMap F F' (toBiratCat P G) biratFT biratDegEq
          (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root)).obj
          (idxOne P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root))) :=
    Under.homMk (biFrHom (biratRtIso F F' X.obj Y.root) hbA
        (biratRtIso F F' Y.obj X.root) hbB hdb)
      (WideSubcategory.hom_ext _ (Prod.ext
        ((Category.id_comp _).trans ((Category.comp_id _).symm.trans
          (congrArg (fun t => biratRtIso F F' X.obj Y.root ≫ t)
            ((toBiratCat P G).map_id _).symm)))
        ((Category.id_comp _).trans ((Category.comp_id _).symm.trans
          (congrArg (fun t => biratRtIso F F' Y.obj X.root ≫ t)
            ((toBiratCat P G).map_id _).symm)))))
  have ht : idxTransport (biratPre P G) F' hv
      (biratRtIso F F' X.obj Y.root ≫ (toBiratCat P G).map χ
        ≫ @inv _ _ _ _ (biratRtIso F F' Y.obj X.root)
          (biratRtIso_isIso F F' Y.obj X.root))
      = (toBiratCat P G).map χ :=
    frobTransport_eq _ _ _ _ _ _ _ (by
      have h1 : ((toBiratCat P G).map χ
            ≫ @inv _ _ _ _ (biratRtIso F F' Y.obj X.root)
              (biratRtIso_isIso F F' Y.obj X.root))
            ≫ biratRtIso F F' Y.obj X.root
          = (toBiratCat P G).map χ :=
        (Category.assoc _ _ _).trans
          ((congrArg (fun t => (toBiratCat P G).map χ ≫ t)
            (IsIso.inv_hom_id (biratRtIso F F' Y.obj X.root))).trans (Category.comp_id _))
      exact (Category.assoc _ _ _).trans
        (congrArg (fun t => biratRtIso F F' X.obj Y.root ≫ t) h1))
  refine Eq.trans (omegaMap_mk F F' (X := X) (Y := Y)
    (idxOne P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root)) χ) ?_
  refine Eq.trans (congrArg (HomPf.mk _) ht.symm) ?_
  exact HomPf.mk_map hv _

variable (F) in
/-- ★共役の計算(圏 `𝒞^birat` の中だけの代数)——
`b₁ ≫ u₁ = r₁⁻¹` と `r_k ≫ b_k = s_k` から `b₁ ≫ (u₁ ≫ s_k) ≫ b_k⁻¹ = r₁⁻¹ ≫ r_k`。

★`Ψ(rtExt A 1)` の可逆性を**要求しない**形にしてある(そこが実装上の要点)。 -/
theorem conj_ext_aux {A R1 S1 Rk Sk : BiratCat P G}
    (r1 : A ⟶ R1) [IsIso r1] (b1 : R1 ⟶ S1) (u1 : S1 ⟶ A)
    (h1 : b1 ≫ u1 = inv r1)
    (rk : A ⟶ Rk) (bk : Rk ⟶ Sk) [IsIso bk] (sk : A ⟶ Sk) (tk : rk ≫ bk = sk) :
    b1 ≫ (u1 ≫ sk) ≫ inv bk = inv r1 ≫ rk := by
  subst tk
  calc b1 ≫ (u1 ≫ rk ≫ bk) ≫ inv bk
      = (b1 ≫ u1) ≫ rk ≫ bk ≫ inv bk := by simp only [Category.assoc]
    _ = inv r1 ≫ rk := by
        rw [h1, IsIso.hom_inv_id, Category.comp_id]

variable (F) in
/-- ★★★★★★**`Ω` は標準の Frobenius 型射 `κ` を `κ` に送る**。

★三角形 `rtExt^birat ≫ β = Ψ(rtExt)` を 2 か所(次数 `1` と `k`)で使うだけ。
★次数 `1` の側は「`β₁ ≫ Ψ(rtOneInv) = (rtExt^birat)⁻¹`」の形で使う ——
こうすると `Ψ(rtExt A 1)` の可逆性を経由せずに済む。 -/
theorem omegaMap_pfKappa (F' : FrobenioidCore (biratPre P G)) (A : C) (k : ℕ+) :
    omegaMap F F' (show HomRoot P F (⟨A, k⟩ : PfRootObj P F) ⟨A, 1⟩ from
        pfKappa (F := F) A k)
      = (show HomRoot (biratPre P G) F'
            (omegaObj F F' (⟨A, k⟩ : PfRootObj P F)) (omegaObj F F' ⟨A, 1⟩) from
          pfKappa (F := F') (biratUp P G A) k) := by
  refine Eq.trans (omegaMap_toHomPf F F' (X := (⟨A, k⟩ : PfRootObj P F)) (Y := ⟨A, 1⟩)
    (kappaRep (P := P) (F := F) A k)) ?_
  refine congrArg (toHomPf (F := F')) ?_
  haveI hr1 : IsIso (rtExt P F A 1) := isIso_rtExt_one P F A
  haveI hq1 : IsIso (rtExt (biratPre P G) F' (biratUp P G A) 1) :=
    isIso_rtExt_one (biratPre P G) F' (biratUp P G A)
  -- ★次数 1 の三角形を `β₁ ≫ Ψ(rtOneInv) = (rtExt^birat)⁻¹` の形にする
  have h1 : biratRtIso F F' A 1
        ≫ (toBiratCat P G).map (rtOneInv (P := P) (F := F) A)
      = @inv _ _ _ _ (rtExt (biratPre P G) F' (biratUp P G A) 1) hq1 := by
    refine (@IsIso.inv_eq_of_hom_inv_id _ _ _ _
      (rtExt (biratPre P G) F' (biratUp P G A) 1) hq1
      (biratRtIso F F' A 1 ≫ (toBiratCat P G).map (rtOneInv (P := P) (F := F) A))
      ?_).symm
    have e2 : rtExt P F A 1 ≫ rtOneInv (P := P) (F := F) A = 𝟙 A := IsIso.hom_inv_id _
    have e3 : (toBiratCat P G).map (rtExt P F A 1)
          ≫ (toBiratCat P G).map (rtOneInv (P := P) (F := F) A)
        = 𝟙 ((toBiratCat P G).obj A) :=
      ((toBiratCat P G).map_comp _ _).symm.trans
        ((congrArg (toBiratCat P G).map e2).trans ((toBiratCat P G).map_id A))
    calc rtExt (biratPre P G) F' (biratUp P G A) 1
            ≫ biratRtIso F F' A 1
              ≫ (toBiratCat P G).map (rtOneInv (P := P) (F := F) A)
        = (rtExt (biratPre P G) F' (biratUp P G A) 1 ≫ biratRtIso F F' A 1)
            ≫ (toBiratCat P G).map (rtOneInv (P := P) (F := F) A) :=
          (Category.assoc _ _ _).symm
      _ = (toBiratCat P G).map (rtExt P F A 1)
            ≫ (toBiratCat P G).map (rtOneInv (P := P) (F := F) A) :=
          congrArg (fun t => t ≫ (toBiratCat P G).map (rtOneInv (P := P) (F := F) A))
            (biratRtIso_tri F F' A 1)
      _ = 𝟙 ((toBiratCat P G).obj A) := e3
  -- ★左辺の中身を `Ψ(rtOneInv) ≫ Ψ(rtExt A k)` に分ける
  have hmid : (toBiratCat P G).map (kappaRep (P := P) (F := F) A k)
      = (toBiratCat P G).map (rtOneInv (P := P) (F := F) A)
        ≫ (toBiratCat P G).map (rtExt P F A k) :=
    (toBiratCat P G).map_comp _ _
  refine Eq.trans (congrArg (fun t => biratRtIso F F' A 1 ≫ t
    ≫ @inv _ _ _ _ (biratRtIso F F' A k) (biratRtIso_isIso F F' A k)) hmid) ?_
  exact conj_ext_aux (rtExt (biratPre P G) F' (biratUp P G A) 1)
    (biratRtIso F F' A 1) ((toBiratCat P G).map (rtOneInv (P := P) (F := F) A)) h1
    (rtExt (biratPre P G) F' (biratUp P G A) k) (biratRtIso F F' A k)
    ((toBiratCat P G).map (rtExt P F A k)) (biratRtIso_tri F F' A k)

/-- ★共役の計算(その 2)—— 真ん中に射を 1 本挟んだ形。 -/
theorem conj_mid_aux {A1 S1 A' B' Rk Sk : BiratCat P G}
    (b1 : A1 ⟶ S1) (u1 : S1 ⟶ A') (w1 : A1 ⟶ A') (h1 : b1 ≫ u1 = w1)
    (m : A' ⟶ B') (rk : B' ⟶ Rk) (bk : Rk ⟶ Sk) [IsIso bk] (sk : B' ⟶ Sk)
    (tk : rk ≫ bk = sk) :
    b1 ≫ (u1 ≫ m ≫ sk) ≫ inv bk = w1 ≫ m ≫ rk := by
  subst tk
  subst h1
  simp

variable (F) in
/-- ★★★★★★**`Ω` は `𝒞 → 𝒞^pf` の射を `𝒞^birat → (𝒞^birat)^pf` の射に送る**。

★`Ω ∘ toRootHom = toRootHom ∘ Ψ` —— 1-可換性の射の部分である。 -/
theorem omegaMap_toRootHom (F' : FrobenioidCore (biratPre P G)) {A B : C} (f : A ⟶ B) :
    omegaMap F F' (show HomRoot P F (⟨A, 1⟩ : PfRootObj P F) ⟨B, 1⟩ from
        toRootHom (F := F) f)
      = (show HomRoot (biratPre P G) F'
            (omegaObj F F' (⟨A, 1⟩ : PfRootObj P F)) (omegaObj F F' ⟨B, 1⟩) from
          toRootHom (F := F') ((toBiratCat P G).map f)) := by
  haveI hrA : IsIso (rtExt P F A 1) := isIso_rtExt_one P F A
  haveI hqA : IsIso (rtExt (biratPre P G) F' (biratUp P G A) 1) :=
    isIso_rtExt_one (biratPre P G) F' (biratUp P G A)
  refine Eq.trans (omegaMap_toHomPf F F' (X := (⟨A, 1⟩ : PfRootObj P F)) (Y := ⟨B, 1⟩)
    (@inv _ _ _ _ (rtExt P F A 1) hrA ≫ f ≫ rtExt P F B 1)) ?_
  refine congrArg (toHomPf (F := F')) ?_
  -- ★次数 1 の三角形(`A` 側)を `β ≫ Ψ(rtExt A 1)⁻¹ = (rtExt^birat)⁻¹` の形に
  have h1 : biratRtIso F F' A 1
        ≫ (toBiratCat P G).map (@inv _ _ _ _ (rtExt P F A 1) hrA)
      = @inv _ _ _ _ (rtExt (biratPre P G) F' (biratUp P G A) 1) hqA := by
    refine (@IsIso.inv_eq_of_hom_inv_id _ _ _ _
      (rtExt (biratPre P G) F' (biratUp P G A) 1) hqA
      (biratRtIso F F' A 1
        ≫ (toBiratCat P G).map (@inv _ _ _ _ (rtExt P F A 1) hrA)) ?_).symm
    have e2 : rtExt P F A 1 ≫ @inv _ _ _ _ (rtExt P F A 1) hrA = 𝟙 A :=
      IsIso.hom_inv_id _
    have e3 : (toBiratCat P G).map (rtExt P F A 1)
          ≫ (toBiratCat P G).map (@inv _ _ _ _ (rtExt P F A 1) hrA)
        = 𝟙 ((toBiratCat P G).obj A) :=
      ((toBiratCat P G).map_comp _ _).symm.trans
        ((congrArg (toBiratCat P G).map e2).trans ((toBiratCat P G).map_id A))
    calc rtExt (biratPre P G) F' (biratUp P G A) 1
            ≫ biratRtIso F F' A 1
              ≫ (toBiratCat P G).map (@inv _ _ _ _ (rtExt P F A 1) hrA)
        = (rtExt (biratPre P G) F' (biratUp P G A) 1 ≫ biratRtIso F F' A 1)
            ≫ (toBiratCat P G).map (@inv _ _ _ _ (rtExt P F A 1) hrA) :=
          (Category.assoc _ _ _).symm
      _ = (toBiratCat P G).map (rtExt P F A 1)
            ≫ (toBiratCat P G).map (@inv _ _ _ _ (rtExt P F A 1) hrA) :=
          congrArg (fun t => t ≫ (toBiratCat P G).map
            (@inv _ _ _ _ (rtExt P F A 1) hrA)) (biratRtIso_tri F F' A 1)
      _ = 𝟙 ((toBiratCat P G).obj A) := e3
  -- ★左辺の中身を 3 本に分ける
  have hmid : (toBiratCat P G).map
        (@inv _ _ _ _ (rtExt P F A 1) hrA ≫ f ≫ rtExt P F B 1)
      = (toBiratCat P G).map (@inv _ _ _ _ (rtExt P F A 1) hrA)
        ≫ (toBiratCat P G).map f ≫ (toBiratCat P G).map (rtExt P F B 1) :=
    ((toBiratCat P G).map_comp _ _).trans
      (congrArg (fun t => (toBiratCat P G).map (@inv _ _ _ _ (rtExt P F A 1) hrA) ≫ t)
        ((toBiratCat P G).map_comp _ _))
  refine Eq.trans (congrArg (fun t => biratRtIso F F' A 1 ≫ t
    ≫ @inv _ _ _ _ (biratRtIso F F' B 1) (biratRtIso_isIso F F' B 1)) hmid) ?_
  exact conj_mid_aux (biratRtIso F F' A 1)
    ((toBiratCat P G).map (@inv _ _ _ _ (rtExt P F A 1) hrA)) _ h1
    ((toBiratCat P G).map f) (rtExt (biratPre P G) F' (biratUp P G B) 1)
    (biratRtIso F F' B 1) ((toBiratCat P G).map (rtExt P F B 1))
    (biratRtIso_tri F F' B 1)


set_option maxHeartbeats 1600000 in
variable (F) in
/-- ★★★★★★**`Ω` は「根を上げた射」`rootMap` を `rootMap` に送る**。

★★`rootMap` は `κ` との合成で**一意に決まる**(`rootMap_ext`、`κ` が mono)。
`Ω` は関手で、`κ` と `toRootHom` を正しく写す(上の 2 本)ので、
**方程式 `rootMap ≫ κ_B = κ_A ≫ [f]` を送るだけ**で済む。 -/
theorem omegaMap_rootMap (F' : FrobenioidCore (biratPre P G))
    (hfi : IsOfFrobeniusIsotropicType P)
    (hfiB : IsOfFrobeniusIsotropicType (biratPre P G))
    {A B : C} (f : A ⟶ B) (k : ℕ+) :
    omegaMap F F' (show HomRoot P F (⟨A, k⟩ : PfRootObj P F) ⟨B, k⟩ from
        rootMap (F := F) hfi f k)
      = (show HomRoot (biratPre P G) F'
            (omegaObj F F' (⟨A, k⟩ : PfRootObj P F)) (omegaObj F F' ⟨B, k⟩) from
          rootMap (F := F') hfiB ((toBiratCat P G).map f) k) := by
  have hR : compRoot (biratPre P G) F'
        (rootMap (F := F') hfiB ((toBiratCat P G).map f) k)
        (pfKappa (F := F') (biratUp P G B) k)
      = compRoot (biratPre P G) F' (pfKappa (F := F') (biratUp P G A) k)
        (toRootHom (F := F') ((toBiratCat P G).map f)) :=
    rootMap_spec (F := F') hfiB ((toBiratCat P G).map f) k
  have s1 := congrArg (compRoot (biratPre P G) F'
    (omegaMap F F' (show HomRoot P F (⟨A, k⟩ : PfRootObj P F) ⟨B, k⟩ from
      rootMap (F := F) hfi f k))) (omegaMap_pfKappa F F' B k).symm
  have s2 := (omegaMap_comp F F'
    (show HomRoot P F (⟨A, k⟩ : PfRootObj P F) ⟨B, k⟩ from rootMap (F := F) hfi f k)
    (show HomRoot P F (⟨B, k⟩ : PfRootObj P F) ⟨B, 1⟩ from pfKappa (F := F) B k)).symm
  have s3 := congrArg (omegaMap F F' (X := (⟨A, k⟩ : PfRootObj P F)) (Y := ⟨B, 1⟩))
    (rootMap_spec (F := F) hfi f k)
  have s4 := omegaMap_comp F F'
    (show HomRoot P F (⟨A, k⟩ : PfRootObj P F) ⟨A, 1⟩ from pfKappa (F := F) A k)
    (show HomRoot P F (⟨A, 1⟩ : PfRootObj P F) ⟨B, 1⟩ from toRootHom (F := F) f)
  have s5 := congrArg (fun t => compRoot (biratPre P G) F' t
      (omegaMap F F' (show HomRoot P F (⟨A, 1⟩ : PfRootObj P F) ⟨B, 1⟩ from
        toRootHom (F := F) f))) (omegaMap_pfKappa F F' A k)
  have s6 := congrArg (compRoot (biratPre P G) F' (pfKappa (F := F') (biratUp P G A) k))
    (omegaMap_toRootHom F F' f)
  exact rootMap_ext (F := F') hfiB _ _
    (((((s1.trans s2).trans s3).trans s4).trans s5).trans (s6.trans hR.symm))

variable (F) in
/-- ★★★★★★★**`biratDescFunctor` の入力 `hΩ`** ——
`Ω` は `𝒞^pf` の co-angular pre-step を `(𝒞^birat)^pf` の**同型**に送る。

★co-angularity は使わない(上の測定)。 -/
theorem omegaMap_isIso_of_preStep (F' : FrobenioidCore (biratPre P G))
    (hiso : ∀ W : C, IsIsotropic P W) {X Y : PfRootObj P F} (f : HomRoot P F X Y)
    (h : IsPreStep (pfRootPre P F) f) :
    IsIso (show omegaObj F F' X ⟶ omegaObj F F' Y from omegaMap F F' f) := by
  have hisoB : ∀ W : BiratCat P G, IsIsotropic (biratPre P G) W :=
    birat_isOfIsotropicType (G := G) hiso
  have hfiB : IsOfFrobeniusIsotropicType (biratPre P G) := fun A =>
    ⟨A, 𝟙 A, isFrobeniusType_of_isIso (biratPre P G) (𝟙 A), hisoB A⟩
  obtain ⟨Z, φ, rfl⟩ := HomPf.exists_rep f
  have hps : IsPreStep P φ := (isPreStep_mk_iff (X := X) (Y := Y) Z φ).mp h
  rw [omegaMap_mk]
  refine pfRoot_isOfIsotropicType (F := F') hfiB _ _ _ ?_ ?_
  · exact (isIsometric_mk_iff _ _).mpr (birat_isIsometric_all G _)
  · refine (isPreStep_mk_iff _ _).mpr ⟨?_, ?_⟩
    · show (biratPre P G).degFr ((toBiratCat P G).map φ) = 1
      rw [birat_degFr_map]
      exact hps.1
    · exact (birat_isBaseIsomorphism_iff (P := P) (G := G) φ).mpr hps.2

variable (F) in
/-- ★★★★★★★**`Ω` は `coaPreProp` を同型に送る**(`biratDescFunctor` の入力の形)。 -/
theorem omegaFunctor_isIso_of_coaPre (F' : FrobenioidCore (biratPre P G))
    (hiso : ∀ W : C, IsIsotropic P W) {X Y : PfRootObj P F} (f : X ⟶ Y)
    (h : coaPreProp (pfRootPre P F) f) :
    IsIso ((omegaFunctor F F').map f) :=
  omegaMap_isIso_of_preStep F F' hiso f h.2

/-! ## ★10. `Θ : (𝒞^pf)^birat ⥤ (𝒞^birat)^pf` -/

variable (F) in
/-- ★★★★★★★**[FrdI] Proposition 5.5, (ii) の関手** ——

```
Θ : (𝒞^pf)^birat ⥤ (𝒞^birat)^pf
```

★★`Prop44Univ.lean` の **`𝒞^birat` の普遍性**(`biratDescFunctor`)に
`Ω` を流すだけ —— 合成の coherence は普遍性の側で済んでいる。 -/
noncomputable def thetaFunctor (F' : FrobenioidCore (biratPre P G))
    (Gpf : Frobenioid (pfRootPre P F)) (hiso : ∀ W : C, IsIsotropic P W) :
    BiratCat (pfRootPre P F) Gpf ⥤ PfRootObj (biratPre P G) F' :=
  biratDescFunctor (G := Gpf) (omegaFunctor F F')
    (fun {_ _} f h => omegaFunctor_isIso_of_coaPre F F' hiso f h)

@[simp] theorem thetaFunctor_obj (F' : FrobenioidCore (biratPre P G))
    (Gpf : Frobenioid (pfRootPre P F)) (hiso : ∀ W : C, IsIsotropic P W)
    (X : BiratCat (pfRootPre P F) Gpf) :
    (thetaFunctor F F' Gpf hiso).obj X
      = omegaObj F F' (biratDown (pfRootPre P F) Gpf X) := rfl

/-- ★★★★★**`Θ` は対象の上で全射** —— `⟨A,n⟩ = Θ ⟨A^birat, n⟩` の形。

★★これで `Proposition 5.5, (ii)` の**本質的全射性は無料**である
(源が `(𝒞^pf)^birat` 全体なので、根 1 に閉じ込められていた
`biratPfFunctor` の難点が消える)。 -/
theorem thetaFunctor_obj_surjective (F' : FrobenioidCore (biratPre P G))
    (Gpf : Frobenioid (pfRootPre P F)) (hiso : ∀ W : C, IsIsotropic P W)
    (Y : PfRootObj (biratPre P G) F') :
    ∃ X : BiratCat (pfRootPre P F) Gpf, (thetaFunctor F F' Gpf hiso).obj X = Y := by
  refine ⟨show BiratCat (pfRootPre P F) Gpf from
    ({ obj := biratDown P G Y.obj, root := Y.root } : PfRootObj P F), ?_⟩
  rfl


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

/-- ★★★★★★locator —— `Proposition 5.5, (ii)` の `Ω` は**関手**である。 -/
def omegaFunctor.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — Ω : 𝒞^pf ⥤ (𝒞^birat)^pf(関手)",
    sectionId := "frdi-prop-5-5" }

/-- ★★★★★★locator —— `Proposition 5.5, (ii)` の `Ω` は
**co-angular pre-step を同型に送る**(`𝒞^birat` の普遍性の入力)。 -/
def omegaFunctor_isIso_of_coaPre.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — Ω は co-angular pre-step を同型に送る",
    sectionId := "frdi-prop-5-5" }

/-- ★★★★★★★locator —— `Proposition 5.5, (ii)` の関手
`Θ : (𝒞^pf)^birat ⥤ (𝒞^birat)^pf` と、その**対象の上の全射性**。 -/
def thetaFunctor.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — Θ : (𝒞^pf)^birat ⥤ (𝒞^birat)^pf(対象は全射)",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
