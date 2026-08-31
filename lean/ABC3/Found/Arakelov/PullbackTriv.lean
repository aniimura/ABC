/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.AMetricHom
import ABC3.Found.Arakelov.PicLTPull
import ABC3.Meta.Claim

/-!
# 引き戻した層の**自明化**を名前にする（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

## ★★★★★★★`Definition 1.1, (ii)` —— 台の側は**すべて在庫**である

2026-08-28 に台帳 `arakelov-pullback-monoidal` を立てたのは**在庫確認の怠り**であった。
★以下はすべて `Found/Arakelov/` に `sorry` 無しで在る:

| 在庫 | 内容 |
|---|---|
| `pullbackPre f` | スキームの射に沿った**前層**加群の引き戻し |
| `pullbackPreOplax` | ★`(pullbackPre f).OplaxMonoidal` |
| `pullbackDelta` / `pullbackEta` | `f^*(P⊗Q) ⟶ f^*P ⊗ f^*Q` と `f^*𝒪 ⟶ 𝒪` |
| `isLocallyTrivial_pullbackPre` | ★★局所自明性は引き戻しで保たれる |
| `bcIso` | ★★★Beck–Chevalley `(f|)^*(L|_W) ≅ (f^*L)|_{f⁻¹W}` |
| `pullbackOnUnitIso` | `(f|)^* 𝟙_ ≅ 𝟙_` |
| `pullbackUnitIso` / `pullbackSheafTensorIso'` | 層の段の単位・テンソル |

## ★★★★★★本ファイルがすること

`isLocallyTrivial_pullbackPre` の証明の中に**自明化の輸送が calc として書かれている**が、
名前が付いていないので**外から使えない**。★本ファイルはそれを `pullTrivOfBase` として
取り出す——計量を運ぶ段（次のブロック）はこれを名前で呼ぶ。

## ★★★★★★★★★次のブロックの設計（明示）

`AMetric` の引き戻しは `LocalMetric.tensor`（`TensorMetric.lean`）と**同じ型**である:

1. ★**チャート**: `V : X.Opens`・`p` に対し、`L|_W ≅ 𝟙_` となる `W : Y.Opens` と
   `V' ≤ V ⊓ f⁻¹W` で `p ⁻¹ᵁ V' = ⊤` なるものを取る
   （存在は `L.triv` ＋ `transportSieve`）。
2. ★★**候補値**:

       `‖evalOn p V' (transUnit (f^*L) V' (pullTrivOfBase …) (e|_{V'}))‖⁻¹`
         `· L.metric.h W eW (p ≫ f)`

3. ★★★**チャート独立性**の核は
   「`transUnit` が輸送 `pullTrivOfBase` と可換」——すなわち
   `transUnit (f^*L) V' (pullTrivOfBase eW₁) (pullTrivOfBase eW₂)`
   が `transUnit L W eW₁ eW₂` の**引き戻し**であること。
   ★これは `LocalMetric.tensor` の `transUnit_tensorTriv` に当たる補題である。
4. ★★★★あとは `chartH_triv_indep` / `chartH_shrink` / `chartH_indep` と
   同じ 3 段を書けばよい。

## ★★★★★★★★★段 3 の**道具は特定できた**（2026-08-28）

`transUnit F V e e' = unitEnd V (e.symm ≪≫ e').hom` は **`rfl`** である
（本ファイルの `transUnit_eq_unitEnd`）。★`unitEnd` は**環準同型**
`End(𝟙_|_V) →+* Γ(X,V)` で、在庫に `unitMul_unitEnd` / `unitEnd_unitMul` がある。

★★したがって段 3 は次に翻訳される:

> `pullTrivOfBase` による共役 `Ψ ↦ B⁻¹ ≫ (f|)^*Ψ ≫ B` が
> `End(𝟙_|_W) → End(𝟙_|_{V'})` の上で**環準同型 `ρ` を実現する**こと。

★★★**その翻訳自体は本ファイルで済んだ**（`pullTrivOfBase_comp`）——
`bcIso` が両側で打ち消し合うので、2 つの輸送した自明化の「比」は
**元の「比」の輸送の共役**である。

★★★★したがって**残っているのは 1 本だけ**である:

> `unitEnd V' ((pullUnitIsoOn).inv ≫ G.map Ψ ≫ (pullUnitIsoOn).hom)`
>   `= ρ (unitEnd W Ψ)`   （`G = pullTransportFun`、`ρ : Γ(Y,W) → Γ(X,V')`）

★★★★★**その相方は本ファイルで済んだ**——`evalOn_pullback`:

    `evalOn p V' (ρ u) = evalOn (p ≫ f) W u`

★つまり「関数の側」は閉じており、残るのは「自明化の側」1 本だけである。

★計算の道具も在庫にある:

| 在庫 | 役割 |
|---|---|
| `unitMul_unitEnd` | `Ψ = unitMul W (unitEnd W Ψ)` ——`Ψ` を単元倍に還元 |
| `isoHomUnitGenOn` | ★**制限した site でも `unit` は生成元を生成元へ送る** |
| `pullbackOnCorepresentableBy` | `(f|)^*` の射への作用を `freeYonedaEquiv` で計算する |
| `unitMul_app_apply` / `unitMul_res` | `unitMul` の値と制限 |

## ★★★★★★★★★★**制限の層も剥がれた**——残りは `V'` を含まない 1 本

`unitEnd_restrictOn`（`unitEnd` は制限と可換）と
`unitEnd_pullTransport_reduce` により、残る等式は

    `pullUnitEnd f W Ψ = f.c.app (op W) (unitEnd W Ψ)`

**ただ 1 本**である（`V'` も `L` も現れない。`rfl` ではない、2026-08-28 実測）。

### ★★★★★★★★その証明の筋（導出済み・未実装）

★随伴で移すのが正しい道である:

1. `pullbackPreOn f W ⊣ pushforward (phiOn f W)` の転置で
   `Hom(G' 𝟙_|_W, B) ≅ Hom(𝟙_|_W, F_* B)`。
2. `Hom(𝟙_|_W, M) ≅ M.obj ⟨W,𝟙⟩`（`freeYonedaTermIso` ＋ `freeYonedaEquiv`）。
   ★`(F_* B).obj ⟨W,𝟙⟩ = B.obj ⟨f⁻¹W,𝟙⟩` は**係数の制限**なので、
   台の集合は同じで、`Γ(Y,W)` 作用が `ρ` を通す。
3. `pullbackOnUnitIso.hom` の転置 `η : 𝟙_|_W ⟶ F_* 𝟙_|_{f⁻¹W}` は
   生成元を生成元へ送る（`isoHomUnitGenOn`）ので `η(1) = 1`。
4. ★★したがって `η(u) = η(u • 1) = ρ(u) • η(1) = ρ u`
   ——**`η` の `Γ(Y,W)`-線型性（＝係数制限）が `ρ` を出す**。

★★★これが「なぜ `ρ` が現れるか」の答えである。★実装は残っている。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory TopologicalSpace MonoidalCategory Opposite

/-- ★★★★★★**引き戻した層の自明化**。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

`L|_W ≅ 𝟙_`（`Y` の上）から `(f^*L)|_{V'} ≅ 𝟙_`（`X` の上、`V' ≤ f⁻¹W`）を作る。

★機構は `isLocallyTrivial_pullbackPre` の証明にある 4 段そのもの:
Beck–Chevalley（`bcIso`）→ `eW` の輸送 → `(f|)^*𝟙_ ≅ 𝟙_`（`pullbackOnUnitIso`）→
制限の推移律。★★本ファイルはそれに**名前を付ける**——
計量を運ぶ段が外から呼べるようにするためである。 -/
noncomputable def pullTrivOfBase {X Y : Scheme.{0}} (f : X ⟶ Y) (L : Y.PresheafOfModules)
    (W : Y.Opens) (eW : (restrictPresheafFunctor Y W).obj L ≅ 𝟙_ (PresheafModulesOn Y W))
    {V' : X.Opens} (hV'W : V' ≤ (Opens.map f.base).obj W) :
    (restrictPresheafFunctor X V').obj ((pullbackPre f).obj L)
      ≅ 𝟙_ (PresheafModulesOn X V') :=
  (Iso.refl _ : (restrictPresheafFunctor X V').obj ((pullbackPre f).obj L)
      ≅ (restrictOnFunctor hV'W).obj
          ((restrictPresheafFunctor X ((Opens.map f.base).obj W)).obj ((pullbackPre f).obj L)))
    ≪≫ (restrictOnFunctor hV'W).mapIso (bcIso f W L).symm
    ≪≫ (restrictOnFunctor hV'W).mapIso ((pullbackPreOn f W).mapIso eW)
    ≪≫ (restrictOnFunctor hV'W).mapIso (pullbackOnUnitIso f W)
    ≪≫ (Iso.refl _ : (restrictOnFunctor hV'W).obj
          (𝟙_ (PresheafModulesOn X ((Opens.map f.base).obj W))) ≅ 𝟙_ (PresheafModulesOn X V'))

/-- ★★★★★★**遷移単元は `unitEnd`（単位対象の自己射環）そのものである**。

★`transUnit F V e e' = unitEnd V (e.symm ≪≫ e').hom` は **`rfl`**（2026-08-28 実測）。

★★これが効く理由: `unitEnd` は**環準同型** `End(𝟙_|_V) →+* Γ(X,V)` であり、
在庫に `unitMul_unitEnd`（`φ = unitMul (unitEnd φ)`）と
`unitEnd_unitMul`（`unitEnd (unitMul c) = c`）がある。
★★★したがって「遷移単元が輸送と可換」は
**「輸送が `End(𝟙_)` の上で環準同型 `ρ` を実現する」**に翻訳される。 -/
theorem transUnit_eq_unitEnd {X : Scheme.{0}} (F : X.PresheafOfModules) (V : X.Opens)
    (e e' : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V)) :
    transUnit F V e e' = unitEnd V (e.symm ≪≫ e').hom := rfl

/-! ## ★★★★★★★輸送を**関手**として見る -/

/-- ★**輸送の関手** `PresheafModulesOn Y W ⥤ PresheafModulesOn X V'`。 -/
noncomputable abbrev pullTransportFun {X Y : Scheme.{0}} (f : X ⟶ Y) (W : Y.Opens)
    {V' : X.Opens} (hV'W : V' ≤ (Opens.map f.base).obj W) :
    PresheafModulesOn Y W ⥤ PresheafModulesOn X V' :=
  pullbackPreOn f W ⋙ restrictOnFunctor hV'W

/-- ★★輸送は単位を単位へ送る。 -/
noncomputable def pullUnitIsoOn {X Y : Scheme.{0}} (f : X ⟶ Y) (W : Y.Opens)
    {V' : X.Opens} (hV'W : V' ≤ (Opens.map f.base).obj W) :
    (pullTransportFun f W hV'W).obj (𝟙_ (PresheafModulesOn Y W))
      ≅ 𝟙_ (PresheafModulesOn X V') :=
  (restrictOnFunctor hV'W).mapIso (pullbackOnUnitIso f W)

/-- ★★★★★★★★**2 つの輸送した自明化の「比」は、元の「比」の輸送の共役である**。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

    `(pullTrivOfBase eW).symm ≪≫ (pullTrivOfBase eW')`
      `= (pullUnitIsoOn).symm ≪≫ G.mapIso (eW.symm ≪≫ eW') ≪≫ pullUnitIsoOn`

★機構は圏の代数だけである——`bcIso` は両側に現れて**打ち消し合う**。
★★これで「遷移単元が輸送と可換」は
**「`unitEnd ∘ 共役 ∘ G.map` が環準同型 `ρ` を実現する」**に落ちる。 -/
theorem pullTrivOfBase_comp {X Y : Scheme.{0}} (f : X ⟶ Y) (L : Y.PresheafOfModules)
    (W : Y.Opens) (eW eW' : (restrictPresheafFunctor Y W).obj L ≅ 𝟙_ (PresheafModulesOn Y W))
    {V' : X.Opens} (hV'W : V' ≤ (Opens.map f.base).obj W) :
    (pullTrivOfBase f L W eW hV'W).symm ≪≫ (pullTrivOfBase f L W eW' hV'W)
      = (pullUnitIsoOn f W hV'W).symm
        ≪≫ (pullTransportFun f W hV'W).mapIso (eW.symm ≪≫ eW')
        ≪≫ pullUnitIsoOn f W hV'W := by
  ext1
  simp [pullTrivOfBase, pullUnitIsoOn, pullTransportFun]
  cat_disch

/-! ## ★★★★★評価は引き戻しと可換 -/

/-- ★`p` が `V' ≤ f⁻¹W` に入るなら `p ≫ f` は `W` に入る。 -/
theorem comp_preimage_eq_top_of_le {X Y : Scheme.{0}} (f : X ⟶ Y)
    {p : Spec (CommRingCat.of ℂ) ⟶ X} {W : Y.Opens} {V' : X.Opens}
    (hV'W : V' ≤ (Opens.map f.base).obj W) (hp : p ⁻¹ᵁ V' = ⊤) :
    (p ≫ f) ⁻¹ᵁ W = ⊤ := by
  have h : (p ≫ f) ⁻¹ᵁ W = p ⁻¹ᵁ ((Opens.map f.base).obj W) := rfl
  rw [h]
  exact preimage_eq_top_of_le hV'W hp

/-- ★★★★★★**評価は関数の引き戻しと可換である**:

    `evalOn p V' (ρ u) = evalOn (p ≫ f) W u`

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

★★これが「引き戻した計量が `p` で見ると元の計量を `p ≫ f` で見たものになる」の
**関数の側**である。★機構は mathlib の `Scheme.Hom.appLE_comp_appLE` と
在庫の `evalOn_restrict` だけである。 -/
theorem evalOn_pullback {X Y : Scheme.{0}} (f : X ⟶ Y)
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (W : Y.Opens)
    {V' : X.Opens} (hV'W : V' ≤ (Opens.map f.base).obj W) (hp : p ⁻¹ᵁ V' = ⊤)
    (hpW : (p ≫ f) ⁻¹ᵁ W = ⊤) (u : (Γ(Y, W) : Type)) :
    evalOn p V' hp (X.presheaf.map (homOfLE hV'W).op (f.c.app (op W) u))
      = evalOn (p ≫ f) W hpW u := by
  rw [evalOn_restrict p hV'W hp]
  show (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom
      ((p.appLE ((Opens.map f.base).obj W) ⊤ _).hom ((f.c.app (op W)).hom u)) = _
  congr 1

/-! ## ★★★★★★★制限の層を剥がす -/

/-- ★★★★**`unitEnd` は制限と可換である**。

★機構は `unitMul_unitEnd` → `unitMul_res` → `unitEnd_unitMul` の 3 段。 -/
theorem unitEnd_restrictOn {X : Scheme.{0}} {U V' : X.Opens} (h : V' ≤ U)
    (Φ : End (𝟙_ (PresheafModulesOn X U))) :
    unitEnd V' ((restrictOnFunctor h).map Φ)
      = X.presheaf.map (homOfLE h).op (unitEnd U Φ) := by
  conv_lhs => rw [← unitMul_unitEnd U Φ]
  rw [unitMul_res h (unitEnd U Φ), unitEnd_unitMul]

/-- ★★**引き戻しが `End(𝟙_)` の上に誘導する写像**（制限を含まない形）。 -/
noncomputable def pullUnitEnd {X Y : Scheme.{0}} (f : X ⟶ Y) (W : Y.Opens)
    (Ψ : End (𝟙_ (PresheafModulesOn Y W))) :
    (Γ(X, (Opens.map f.base).obj W) : Type) :=
  unitEnd ((Opens.map f.base).obj W)
    ((pullbackOnUnitIso f W).inv ≫ (pullbackPreOn f W).map Ψ ≫ (pullbackOnUnitIso f W).hom)

/-- ★★★★★★**輸送の共役から制限の層が剥がれる**。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

★これで「遷移単元が輸送と可換」の残りは **`V'` を含まない 1 本**になった:

    `pullUnitEnd f W Ψ = f.c.app (op W) (unitEnd W Ψ)`

★★機構は `restrictOnFunctor` の関手性と `unitEnd_restrictOn` だけである。 -/
theorem unitEnd_pullTransport_reduce {X Y : Scheme.{0}} (f : X ⟶ Y) (W : Y.Opens)
    {V' : X.Opens} (hV'W : V' ≤ (Opens.map f.base).obj W)
    (Ψ : End (𝟙_ (PresheafModulesOn Y W))) :
    unitEnd V' ((pullUnitIsoOn f W hV'W).inv ≫ (pullTransportFun f W hV'W).map Ψ
        ≫ (pullUnitIsoOn f W hV'W).hom)
      = X.presheaf.map (homOfLE hV'W).op (pullUnitEnd f W Ψ) := by
  have hcomp : (pullUnitIsoOn f W hV'W).inv ≫ (pullTransportFun f W hV'W).map Ψ
      ≫ (pullUnitIsoOn f W hV'W).hom
      = (restrictOnFunctor hV'W).map ((pullbackOnUnitIso f W).inv
          ≫ (pullbackPreOn f W).map Ψ ≫ (pullbackOnUnitIso f W).hom) := by
    simp only [pullUnitIsoOn, pullTransportFun, Functor.mapIso_inv, Functor.mapIso_hom,
      Functor.comp_map, Functor.map_comp]
    rfl
  rw [hcomp, unitEnd_restrictOn]
  rfl

/-! ### ★出典の紐付け(`.src`) -/

def pullTrivOfBase.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(引き戻した層の自明化——輸送を名前にしたもの)",
    sectionId := "genell-def-1-1-ii" }

def pullTrivOfBase_comp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(輸送した 2 つの自明化の比は、元の比の輸送の共役であること)",
    sectionId := "genell-def-1-1-ii" }

def unitEnd_restrictOn.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(unitEnd は制限と可換であること)",
    sectionId := "genell-def-1-1-i" }

def unitEnd_pullTransport_reduce.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(輸送の共役から制限の層が剥がれること)",
    sectionId := "genell-def-1-1-ii" }

def evalOn_pullback.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(評価は関数の引き戻しと可換であること)",
    sectionId := "genell-def-1-1-ii" }

def transUnit_eq_unitEnd.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(遷移単元は単位対象の自己射環 unitEnd そのものであること)",
    sectionId := "genell-def-1-1-i" }

def pullTrivOfBase.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "bcIso(Beck–Chevalley)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.bcIso") 3,
    .citation "[ABC3]" "pullbackOnUnitIso((f|)^* 𝟙_ ≅ 𝟙_)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.pullbackOnUnitIso") 3,
    .citation "[ABC3]" "isLocallyTrivial_pullbackPre(局所自明性は引き戻しで保たれる)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.isLocallyTrivial_pullbackPre") 3,
    .implicitStep
      ("★★★★★**訂正(2026-08-28)**: 台帳 arakelov-pullback-monoidal を立てたのは" ++
       "在庫確認の怠りであった。mathlib に Monoidal インスタンスが無いのは事実だが、" ++
       "本プロジェクトは pullbackPreOplax を既に建てており、" ++
       "自明化の輸送も isLocallyTrivial_pullbackPre の証明の中にある") 3,
    .implicitStep
      ("★★次のブロックの設計: AMetric の引き戻しは LocalMetric.tensor と同じ型である。" ++
       "(1) チャート(L|_W ≅ 𝟙_ となる W と V' ≤ V ⊓ f⁻¹W)、" ++
       "(2) 候補値 ‖evalOn(transUnit (f^*L) V' (pullTrivOfBase …) (e|_{V'}))‖⁻¹ · h W eW (p ≫ f)、" ++
       "(3) チャート独立性の核は「transUnit が輸送と可換」" ++
       "(LocalMetric.tensor の transUnit_tensorTriv に当たる補題)。" ++
       "★見積もり 2〜4 ブロック") 3 ]

end ABC3.Found.Arakelov
