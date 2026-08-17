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

★本ファイルは**そこまでは主張しない**。取ったのは第 1 段だけである。

### ★★★器具の記録 —— 易しい側の包含すら組めなかった

`Ideal.map (f.appLE V U h) (I.ideal V) ≤ (I.comap f).ideal U`(易しい側)を
10 回試して**組み上がらなかった**。★数学は 3 行である:
`pullback.condition` を `V` で `app` に落とし、
`s ∈ I.ideal V = ker (I.subschemeι.app V)` を使い、`appLE_comp_appLE` で `U` へ制限する。

★★**摩擦は `CommRingCat.Hom.hom` と `ConcreteCategory.hom` の正規形の争い**である。
`simp only` を当てるたびに**どちらかへ揺れ戻り**、`rw [h0]` が当たらない。
`CommRingCat.comp_apply` / `ConcreteCategory.comp_apply` / `RingHom.comp_apply` /
`CommRingCat.hom_comp` のどれを入れても、別の形に正規化されて次の `rw` が外れた。

★★★**同じ証明が、import の重さで通ったり通らなかったりした**——
軽い import の probe では通り、同じ内容を本ファイルに置くと通らない。
★これは本日 5 例目の「同じ数学が書き方で通ったり通らなかったりする」であり、
**今までで最も悪質な型**である(書き方ではなく**環境**で変わる)。

★次に取るときは `ConcreteCategory` の API に統一して書くこと。
**半端に `sorry` を置かず、取れた第 1 段だけを残す。**
-/

/-! ## ★出典の紐付け(`.src`) -/

def ideal_comap_eq_ker.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (i)(comap のアフィン局所記述——核として書ける所まで)",
    sectionId := "genell-prop-1-4" }

end ABC3.Found.GenEll
