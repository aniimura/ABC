import ABC3.Meta.Claim
import Mathlib.AlgebraicGeometry.IdealSheaf.Functorial

/-!
# [GenEll] Proposition 1.4, (i) の残り —— `comap` の**アフィン局所記述**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

## ★★★`Proposition 1.4, (i)` に残った 1 本を正面から取りにいく

`HeightConstruction.lean` は `ht` の加法性を **`PullbackMul` 1 本**に落とした。
その正体は `Scheme.IdealSheafData.comap_mul`(引き戻しが積を保つこと)であり、
★mathlib には**無い**(2026-08-17 実測。あるのは `comap_sup`——左随伴だから上限は保つ)。

## ★★なぜ自明でないのか

mathlib の定義は

    `comap I f := (pullback.fst f I.subschemeι).ker`

すなわち**ファイバー積の射影の核**である。
★積との両立を見るには**アフィン局所の記述**が要るが、
mathlib にあるのは `ideal_comap_of_isOpenImmersion`(開埋め込みの場合)だけである。

## ★★★本ファイルが取る段

`comap I f = map ⊥ (pullback.fst f I.subschemeι)`(`map_bot` から)であり、
★★`I.subschemeι` は**閉埋め込み**なのでアフィン射、
したがって `pullback.fst f I.subschemeι` も**アフィン射**である(底変換で保たれる)。
★`ideal_map_of_isAffineHom` が使えて、**核として書ける**:

    `(I.comap f).ideal U = RingHom.ker ((pullback.fst f I.subschemeι).app U)`

★★★**これが `comap` のアフィン局所記述の第 1 段**である。
残るのは「その核が `Ideal.map` に等しい」ことで、
そこには `Γ(ファイバー積) = テンソル積`(`pullbackSpecIso`)が要る。

## ★実測(2026-08-17)

- `IsClosedImmersion I.subschemeι` —— instance で出る
- `IsAffineHom (pullback.fst f I.subschemeι)` —— instance で出る
- `QuasiCompact (pullback.fst f I.subschemeι)` —— instance で出る
- `(⊥ : X.IdealSheafData).ideal U = ⊥` —— `rfl`
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory CategoryTheory.Limits

variable {X Y : Scheme}

/-! ## ★★`comap` は `⊥` の押し出しである -/

/-- ★**`comap I f` はファイバー積の射影の核である**(定義の言い換え)。

★`map_bot` を使って `map`(押し出し)の言葉に直す——
そうすると `ideal_map_of_isAffineHom` が使える。 -/
theorem comap_eq_map_bot (I : Y.IdealSheafData) (f : X ⟶ Y) :
    I.comap f = Scheme.IdealSheafData.map ⊥ (pullback.fst f I.subschemeι) := by
  rw [Scheme.IdealSheafData.map_bot]
  rfl

/-! ## ★★★アフィン局所記述(第 1 段) -/

/-- ★★★**`comap` のアフィン開集合上の値は核である**。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

★★機構は 2 つ:
`I.subschemeι` が**閉埋め込み**(⟹ アフィン射)であり、
アフィン射は**底変換で保たれる**ので `pullback.fst f I.subschemeι` もアフィン射。
★したがって `ideal_map_of_isAffineHom` が使える。

★★★**これで `comap` が「不透明なファイバー積」から
「環準同型の核」に変わった。** -/
theorem ideal_comap_eq_ker (I : Y.IdealSheafData) (f : X ⟶ Y) (U : X.affineOpens) :
    (I.comap f).ideal U
      = RingHom.ker ((pullback.fst f I.subschemeι).app U.1).hom := by
  rw [comap_eq_map_bot, Scheme.IdealSheafData.ideal_map_of_isAffineHom]
  show Ideal.comap _ (⊥ : Ideal _) = _
  rw [← RingHom.ker_eq_comap_bot]

/-- ★**空因子の引き戻しでは核は自明**(負の対照)。

★`comap ⊤ f = ⊤` なので核は `⊤` になる——
`ideal_comap_eq_ker` が向きを取り違えていれば、ここが `⊥` になる。 -/
theorem ideal_comap_top_eq_ker (f : X ⟶ Y) (U : X.affineOpens) :
    RingHom.ker ((pullback.fst f (⊤ : Y.IdealSheafData).subschemeι).app U.1).hom = ⊤ := by
  rw [← ideal_comap_eq_ker, Scheme.IdealSheafData.comap_top]
  rfl

/-! ## ★★★易しい側の包含 —— `Ideal.map ≤ comap` -/

/-- ★★**因子の切断はファイバー積の上で消える**。

★機構は `pullback.condition` を `V` で `app` に落とすことだけ。
`s ∈ I.ideal V = ker (I.subschemeι.app V)`(`ker_subschemeι_app`)を使う。

★★★**器具の要点: 射のレベルで分解し、元に降りるのは 1 度だけにする。**
途中で `simp` を当てると `CommRingCat.Hom.hom` と `ConcreteCategory.hom` の間で
正規形が揺れ戻り、次の `rw` が外れる(下の docstring に記録)。 -/
theorem app_pullback_comp_eq_zero (I : Y.IdealSheafData) (f : X ⟶ Y)
    (V : Y.affineOpens) (s : Γ(Y, V.1)) (hs : s ∈ I.ideal V) :
    ((pullback.fst f I.subschemeι ≫ f).app V.1).hom s = 0 := by
  have hcond := Scheme.Hom.congr_app
    (pullback.condition (f := f) (g := I.subschemeι)) V.1
  have hpre : (pullback.fst f I.subschemeι ≫ f) ⁻¹ᵁ V.1
      = pullback.snd f I.subschemeι ⁻¹ᵁ I.subschemeι ⁻¹ᵁ V.1 := by
    rw [pullback.condition]; rfl
  -- ★射のレベルで分解する(元に降りない)
  have hA : (pullback.fst f I.subschemeι ≫ f).app V.1
      = I.subschemeι.app V.1 ≫
        ((pullback.snd f I.subschemeι).app (I.subschemeι ⁻¹ᵁ V.1) ≫
          (pullback f I.subschemeι).presheaf.map (eqToHom hpre).op) := by
    rw [hcond, Scheme.Hom.comp_app]
    exact Category.assoc _ _ _
  -- ★ここで 1 度だけ元に降りる
  rw [hA, CommRingCat.hom_comp, RingHom.comp_apply]
  have h0 : (I.subschemeι.app V.1).hom s = 0 := by
    rw [← RingHom.mem_ker, Scheme.IdealSheafData.ker_subschemeι_app]
    exact hs
  rw [h0, map_zero]

/-- ★★★**引き戻したイデアルは `comap` に含まれる**。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

★★これが `comap` のアフィン局所記述の**易しい側**である:

    `Ideal.map (f.appLE V U h) (I.ideal V) ≤ (I.comap f).ideal U`

★逆向き(`comap ≤ Ideal.map`)にはファイバー積の `Γ` の計算が要る。

★★機構は `appLE_comp_appLE` で `V → U` の制限を合成に直し、
`app_pullback_comp_eq_zero` を当てるだけである。 -/
theorem map_appLE_le_ideal_comap (I : Y.IdealSheafData) (f : X ⟶ Y)
    (U : X.affineOpens) (V : Y.affineOpens) (h : U.1 ≤ f ⁻¹ᵁ V.1) :
    (I.ideal V).map (f.appLE V.1 U.1 h).hom ≤ (I.comap f).ideal U := by
  rw [ideal_comap_eq_ker, Ideal.map_le_iff_le_comap]
  intro s hs
  simp only [Ideal.mem_comap, RingHom.mem_ker]
  -- ★射のレベルで: `appLE f ≫ g.app U = (g ≫ f).appLE V (g⁻¹U)`
  have key : (Scheme.Hom.appLE f V.1 U.1 h) ≫ ((pullback.fst f I.subschemeι).app U.1)
      = (pullback.fst f I.subschemeι ≫ f).appLE V.1
          ((pullback.fst f I.subschemeι) ⁻¹ᵁ U.1) (fun x hx => h hx) := by
    rw [Scheme.Hom.app_eq_appLE, Scheme.Hom.appLE_comp_appLE]
  -- ★元に降りるのは 1 度だけ
  have hdesc : ((pullback.fst f I.subschemeι).app U.1).hom ((f.appLE V.1 U.1 h).hom s)
      = ((pullback.fst f I.subschemeι ≫ f).appLE V.1
          ((pullback.fst f I.subschemeι) ⁻¹ᵁ U.1) (fun x hx => h hx)).hom s := by
    rw [← key, CommRingCat.hom_comp, RingHom.comp_apply]
  rw [hdesc, Scheme.Hom.appLE, CommRingCat.hom_comp, RingHom.comp_apply,
    app_pullback_comp_eq_zero I f V s hs, map_zero]

/-! ## ★★残る段を明示する —— 部品名まで特定した(2026-08-17 夜)

★★★**核が `Ideal.map` に等しい**ことが残る:

    `RingHom.ker ((pullback.fst f I.subschemeι).app U) = (I.ideal V).map (f.appLE V U _)`

### ★実測: 使える部品は mathlib に揃っている

| 部品 | 場所 | 内容 |
|---|---|---|
| `ker_subschemeι_app` | `IdealSheaf/Subscheme.lean` | `ker (I.subschemeι.app U) = I.ideal U` |
| `subschemeι_app_surjective` | 同上 | `I.subschemeι.app U` は**全射** |
| `subschemeObjIso` | 同上 | `Γ(𝒪ₓ/I, U) ≅ Γ(X,U)/I(U)` |
| `pullbackSpecIso` | `AlgebraicGeometry/` | `Γ(ファイバー積) = テンソル積` |
| `comap_comp` | `IdealSheaf/Functorial.lean` | アフィンへの帰着に使う |
| `ideal_comap_of_isOpenImmersion` | 同上 | 開埋め込みへの制限 |

### ★★★筋道(次のセッションはここから)

1. `comap_comp` + `ideal_comap_of_isOpenImmersion` で
   **`X`, `Y` がアフィンな場合に帰着する**
   (`U ↪ X → Y` を `U → V ↪ Y` と分解する)
2. アフィンでは `f = Spec.map φ`、`I.subscheme ≅ Spec (R/I₀)` なので
   `pullbackSpecIso` により `Γ(ファイバー積) ≅ A ⊗_R R/I₀ ≅ A/I₀A`
3. `ker (A → A/I₀A) = I₀·A = Ideal.map φ I₀`
4. `Ideal.map_mul` から **`comap_mul`**

★★これが入れば `Proposition 1.4, (i)` が構成から**無条件に**出る。

★本ファイルが取ったのは**第 1 段(核として書ける)と易しい側の包含**である。

### ★★★器具の記録 —— 「降り方」で通った

★**一度は「10 回試して組めなかった」と書いた。手を変えたら通った。**
数学は 3 行である——`pullback.condition` を `V` で `app` に落とし、
`s ∈ I.ideal V = ker (I.subschemeι.app V)` を使い、`appLE_comp_appLE` で `U` へ制限する。

**通らない書き方**: 途中で `simp only` を当てる。
★`CommRingCat.Hom.hom` と `ConcreteCategory.hom` の間で**正規形が揺れ戻り**、
次の `rw` が外れる。`CommRingCat.comp_apply` / `ConcreteCategory.comp_apply` /
`RingHom.comp_apply` / `CommRingCat.hom_comp` のどれを入れても別の形に正規化される。
★★同じ証明が **import の重さで**通ったり通らなかったりした。

★★★**通る書き方: 射のレベルで分解を完成させ、元に降りるのは 1 度だけにする。**
`Category.assoc` を `rw` ではなく `exact Category.assoc _ _ _` で当て、
最後に `CommRingCat.hom_comp` + `RingHom.comp_apply` で 1 度だけ降りる。
★★**`simp` を 1 度も使わない。**

★これが本日 5 例目の「同じ数学が書き方で通ったり通らなかったりする」であり、
**規則が 1 つ出た**——**圏の言葉で終わらせ、元には最後に 1 度だけ降りる。**

★次に取るときは `ConcreteCategory` の API に統一して書くこと。
**半端に `sorry` を置かず、取れた第 1 段だけを残す。**
-/

/-! ## ★出典の紐付け(`.src`) -/

def ideal_comap_eq_ker.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (i)(comap のアフィン局所記述——核として書ける所まで)",
    sectionId := "genell-prop-1-4" }

end ABC3.Found.GenEll
