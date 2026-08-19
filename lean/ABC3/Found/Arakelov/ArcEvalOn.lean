import ABC3.Found.Arakelov.ArcOverBridge

/-!
# Arakelov (C3) 第 256 ブロック —— ★★★★★**開集合上の切断を評価する**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★§9-297 の障害を**迂回する**

自明化は開集合 `V` の上でしか取れないので、これまでは
`Scheme.Modules.restrict` へ橋を架けようとして
**加群構造が 2 つある**という壁に当たった(§9-297)。

★★★迂回路: **随伴の単位射を `⊤` ではなく `V` で評価する**。

    Γ(V, F) --[η_F の V 成分]--> Γ(p⁻¹V, p^* F)

★`p` が `V` を通れば `p⁻¹V = ⊤` なので、右辺は **`arcFiber p F` そのもの**である。

★★★★**制限関手も引き戻しとの比較も要らない**——`X` の上に留まったままでよい。

## ★★これで何が変わるか

| 旧路(§9-294–297) | 新路(本ブロック) |
|---|---|
| `restrict F V.ι ≅ 𝒪_V` を作る | ★不要 |
| `V.Opens ≌ Over V` を作る | ★不要 |
| 加群構造の二重路 | ★★**現れない** |

★★★自明化 `e` から生成切断 `g ∈ Γ(V, F)` を取り、
`arcEvalOnTop p V h g` がファイバーを張ることを使えばノルムが書ける。

| 定義・定理 | 内容 |
|---|---|
| `arcEvalOn` | ★開集合 `U` 上の切断を `p⁻¹U` の切断へ |
| `preimage_eq_top` | ★`p` が `U` を通れば `p⁻¹U = ⊤` |
| `arcEvalOnTop` | ★★★**`U` 上の切断からファイバーの元へ** |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite

variable {X : Scheme.{0}} (F : X.Modules) (p : Spec (CommRingCat.of ℂ) ⟶ X)

/-- ★**開集合 `U` 上の切断を `p⁻¹U` の切断へ送る**(随伴の単位)。 -/
noncomputable def arcEvalOn (U : X.Opens) (s : (F.val.obj (op U) : Type)) :
    (((Scheme.Modules.pullback p).obj F).val.obj (op (p ⁻¹ᵁ U)) : Type) :=
  (((Scheme.Modules.pullbackPushforwardAdjunction p).unit.app F).val.app (op U)).hom s

/-- ★**`p` が `U` を通るとき `p⁻¹U = ⊤`**。 -/
theorem preimage_eq_top (U : X.Opens) (h : ∀ z, p.base z ∈ U) : p ⁻¹ᵁ U = ⊤ := by
  ext z
  simp only [TopologicalSpace.Opens.coe_top, Set.mem_univ, iff_true]
  exact h z

/-- ★★★**`p` が `U` を通るとき、`U` 上の切断はファイバーの元を定める**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが §9-297 の障害(加群構造の二重路)を迂回する道具である
——`restrict` を経由しないので、構造は 1 つしか現れない。 -/
noncomputable def arcEvalOnTop (U : X.Opens) (h : p ⁻¹ᵁ U = ⊤)
    (s : (F.val.obj (op U) : Type)) : ↥(arcFiber p F) :=
  ((((Scheme.Modules.pullback p).obj F).val.map
    (homOfLE (le_of_eq h.symm)).op)).hom (arcEvalOn F p U s)

/-! ## ★出典の紐付け(`.src`) -/

def arcEvalOnTop.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——開集合上の切断を評価する)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
