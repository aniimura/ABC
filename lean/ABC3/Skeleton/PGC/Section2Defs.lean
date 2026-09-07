import ABC3.Skeleton.PGC.Setup
import ABC3.Interface.PGC.LocalFieldData

/-!
# [pGC] §2 — Proposition 2.1 / 2.2 が語る対象の**定義**

主張の本体は `ABC3/Skeleton/PGC/Section2.lean` にある。
本ファイルはその**定義だけ**を持つ（`Section1Defs.lean` と同じ作法）。

## ★なぜ定義と主張を分けたか（2026-09-07）

`Found/PGC/AbsClosureModules.lean` が `RecoverableAsAddModule` を**使う**。
定義が定理ファイル `Section2.lean` にあると、`Found` 側が
`Section2 → Section1` を引き込むため、
★**Proposition 1.1 の配線（`Section1` → `Found/PGC/ReciprocityDatumIndependence`）で
import 循環になる**：
`Section1 → RDI → AbsClosureModules → Section2 → Section1`。

★名前空間は同じ `ABC3.Skeleton.PGC` なので完全修飾名は 1 文字も変わらない。
-/

namespace ABC3.Skeleton.PGC

open ABC3.Meta ABC3.Interface.PGC

variable {p : ℕ} [Fact p.Prime]

/-! ## Γ_K-加群としての回復可能性(§1 の `RecoverableFromAbsGal` とは別の形)

§1 の3命題(Prop 1.1・1.2・Cor 1.3)はいずれも「単一の値・部分群が回復できる」形
(`AssociatedObject.RecoverableFromAbsGal`)。Proposition 2.1 はこれと異なり、
**加法群としての K̄ の Γ_K-作用込みの構造そのもの**が回復できる、という主張——
「対応する値」ではなく「対応する加法的同型」の存在を述べる。ゆえに新しい形を導入する。 -/

/-- **Γ_K-加群としての回復可能性**。`Obj K` が Γ_K-作用込みの加法群として与えられているとき、
任意の同型 α : Γ_K ≅ Γ_K′ に対し、α と両立する加法的同型 `Obj K ≃+ Obj K′` が存在する
という主張。

## ★★★2026-09-05: 作用のクラスを `SMul` から `DistribMulAction` に直した

原文は「Γ_K-**加群**」と言うが、旧形は作用を **`SMul`**(`one_smul` も
`mul_smul` も `smul_add` も要求しない、公理ゼロのクラス)で受け取っていた。
その形では `prop_2_2` は**偽**である——`Γ_{ℚ_p}` の非可換性
(`Found/PGC/QpNonAbelian.lean`)を使って病的な「作用」
`g • n := if g = g₀ then 0 else n` を作れば反例になる
(`Check/PGC/Prop22Degenerate.lean::prop_2_2_statement_false`、`sorry` 無し)。

`DistribMulAction`(加法的自己同型として作用する)に強めると、この病的な
作用は `MulAction` ですらないので塞がる。`Prop 2.1` が使う `K.closure` への
自然な作用はもちろん `DistribMulAction` を満たす。 -/
def RecoverableAsAddModule (Obj : PAdicLocalField p → Type*)
    [∀ K, AddCommGroup (Obj K)] [∀ K, DistribMulAction K.absGal (Obj K)] : Prop :=
  ∀ {K K' : PAdicLocalField p} (α : ContinuousMulEquiv K.absGal K'.absGal),
    ∃ φ : Obj K ≃+ Obj K', ∀ (g : K.absGal) (x : Obj K),
      φ (g • x) = (α.toMulEquiv g) • (φ x)

/-- §2 冒頭(Proposition 2.1 の準備)で導入する我々自身の定義——原典の項目そのものではない
(bare な `"Proposition 2.1"` と紛れないよう、item 名に注記を付ける)。 -/
def RecoverableAsAddModule.src : Source :=
  { paper := "pGC", pdfPage := 4, item := "Proposition 2.1 (RecoverableAsAddModule)",
    sectionId := "prop-2-1" }

/-- `K.closure`(= K̄)への Γ_K の自然な作用(AlgEquiv としての適用)。

★2026-09-05: `SMul`(公理ゼロ)から `DistribMulAction`(=原文の「Γ_K-加群」)に
強めた。理由は `Check/PGC/Prop22Degenerate.lean` を参照——`SMul` のままだと
`prop_2_2` が**偽**になる。 -/
noncomputable instance closureDistribMulAction (K : PAdicLocalField p) :
    DistribMulAction K.absGal K.closure where
  smul g x := g x
  one_smul _ := rfl
  mul_smul _ _ _ := rfl
  smul_zero g := map_zero g
  smul_add g x y := map_add g x y

end ABC3.Skeleton.PGC
