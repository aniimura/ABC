import ABC3.Found.PGC.LocalTateDualityRoute
import Mathlib.RepresentationTheory.Homological.GroupCohomology.LowDegree
import Mathlib.RepresentationTheory.Homological.GroupCohomology.Functoriality

/-!
# [pGC] Proposition 1.1 —— 双対性の道の「離散側の欠落」を埋める

`Found/PGC/LocalTateDualityRoute.lean` は Proposition 1.1 を
`LocalTateDualityData`(3 つのフィールドからなる仮定)に帰着させた。
そのファイルの測定節は、**離散側に 3 つの欠落**があると名指しした:

* **(a)** `Finite (groupCohomology.H2 A)` が `[Finite G]` 付きでも `failed to synthesize`
* **(b)** `e : G ≃* H` に沿った `H2 (Rep.res e.toMonoidHom A) ≅ H2 A` が束ねられていない
* **(c)** `H²` の inflation-restriction と有限商の塔の colimit が無い

本ファイルは **(a) と (b) を埋め**、**(c) を測り直した**。

## ★本ファイルが無条件に証明したこと(仮定ゼロ)

1. **(a) は埋まった。しかも全次数で。**
   `finiteGroupCohomology : [Finite G] → [Finite ↑A] → Finite ↑(groupCohomology A n)`。
   ★`H2` は `groupCohomology A 2` の `abbrev` なので `Finite ↑(H2 A)` は `inferInstance` になる。
   ★おまけに濃度の上界 `natCard_groupCohomology_le` も付いた
   (`|H^n(G,A)| ≤ |A|^(|G|^n)`)。
2. **(b) は埋まった。しかも全次数で。**
   `groupCohomologyMulEquivIso (e : G ≃* H) (A : Rep k H) (n : ℕ) :
      groupCohomology A n ≅ groupCohomology (Rep.res e.toMonoidHom A) n`。
   ★原典が要求するのは `n = 2` だけだが、証明は次数に一切触れないので全次数で通る。
3. **★★原典の「群論性」は定理である。**
   `discreteCardH2_isGroupTheoretic` —— `LocalTateDualityData.isGroupTheoretic` と
   **同じ形の命題**を、`cardH2` を離散群コホモロジーで定義した場合について証明した。
   ★ただしこれは `LocalTateDualityData` を作るものでは**ない**(下の留保を見よ)。
4. **`charRep`** —— 指標 `ψ : G →* kˣ` から階数 1 の表現 `Rep k G` を作る。
   これが原典の `M`(`Zp`-加群として `Z/pnZ`、`Γ_K` は指標で作用)にあたる。
5. **★(c) の到達点を宣言として残した** —— `inflation` / `inflationH2`。
   ★**次数 2 でも inflation 写像は引ける**(2026-09-04 の記録の訂正)。
   引けないのは**その完全性**と**塔の colimit** である。

## ★★仮定に置いたもの(名指し)

**本ファイルは仮定を 1 つも置いていない。** `structure` も `axiom` も `sorry` も無い。
`LocalTateDualityData` は **作っていない**(作れない理由は下の (c) の測定)。

★**逆に、本ファイルが `LocalTateDualityData` を作らない理由を明示しておく**:
`discreteCardH2` を `cardH2` に代入すると `isGroupTheoretic` は定理として埋まるが、
`cardH2_eq_natCard`(局所 Tate 双対性)は**偽になる**。
副有限群 `Γ_K` の離散(=連続性を課さない)群コホモロジーは双対性を満たさないからである。
★したがって「`LocalTateDualityData.ofDiscrete` を作る」ことは**しない**。
それは偽の仮説から作る空虚な構成になる。

## ★★★(c) の測定(2026-09-07 実施)

★測り方は 2 つ。索引 grep と REPL の `#check` である。

```
awk -F'\t' '$2 ~ /^groupCohomology\.(H2|cocycles|map|congr)/ {print $2"\t"$4}' \
  .cache/mathlib-index.txt
awk -F'\t' '$3 ~ /RepresentationTheory\/Continuous\// {print $2"\t"$3}' .cache/mathlib-index.txt
awk -F'\t' '$3 ~ /ContinuousCohomology/ {print $2"\t"$3"\t"$4}' .cache/mathlib-index.txt
```

REPL(`ABC3.Found.PGC.LocalTateDualityRoute` +
`Mathlib.RepresentationTheory.Homological.GroupCohomology.{LowDegree,Functoriality,LongExactSequence}`
+ `Mathlib.RepresentationTheory.Continuous.Basic`
+ `Mathlib.Algebra.Category.ContinuousCohomology.Basic`)で `#check` した結果:

| 名前 | 結果 |
|---|---|
| `groupCohomology.infNatTrans` | ★**在る。しかも全次数**。`(S : Subgroup G) [S.Normal] (n : ℕ) : quotientToInvariantsFunctor k S ⋙ functor k (G ⧸ S) n ⟶ functor k G n` |
| `groupCohomology.H1InfRes` / `H1InfRes_exact` | 在る(★**次数 1 だけ**) |
| `groupCohomology.H2InfRes` | ★`Unknown constant` |
| `groupCohomology.infNatTrans_exact` | ★`Unknown constant` |
| `groupCohomology.colimitIso` / `infRes` | ★`Unknown constant` |
| `continuousCohomology` | 在る。`Action (TopModuleCat R) G ⥤ TopModuleCat R` |
| `ContinuousCohomology.continuousCohomologyZeroIso` | 在る(★**次数 0 だけ**。`≅ invariants`) |
| `ContinuousCohomology.H2` / `continuousCohomology.colimitIso` | ★`Unknown constant` |
| `ContRepresentation` | 在る(`RepresentationTheory/Continuous/Basic.lean`, 57 宣言) |
| `ContRepresentation.cohomology` | ★`Unknown constant` |

★★**結論(2026-09-07 の再測定)**: **(c) は今の mathlib でも届かない。** ただし
2026-09-04 の記録は **1 点だけ古かった**ので訂正する:

* ★**訂正**: 「inflation は次数 1 だけ」は不正確だった。**inflation 写像 `infNatTrans` は
  全次数で在る**。無いのはその**完全性/同型性**(次数 2 以上)と**塔の colimit** である。
* ★**`ContRepresentation`(57 宣言)は cohomology を 1 つも持たない。**
  中身は `ContIntertwiningMap` / `Equiv` / `coind` / `restrict` / `trivial` ——
  すなわち**連続表現の圏の対象と射だけ**である。「連続コホモロジーの土台」であって
  「連続コホモロジー」ではない。
* ★**`continuousCohomology` は `ContRepresentation` を使っていない。**
  `Action (TopModuleCat R) G` の上に立っており、両者を繋ぐ橋は無い。
  同定されているのは次数 0(`invariants`)だけである。
* ★**`LongExactSequence`(19 宣言)は同じ群 `G` の短完全列に対する `δ` である。**
  inflation とは方向が違うので (c) には効かない。

★**(c) を埋めるために名指しで足りないもの(4 つ)**:
1. `H²` の inflation-restriction 完全列(`H1(G/S,A^S) → H1(G,A) → H1(S,A)^{G/S} → H2(G/S,A^S) → H2(G,A)`)
2. 開正規部分群の有向系についての `colim_S H^n(G/S, A^S) ≅ H^n_cont(G, A)`
3. `ContRepresentation` と `Action (TopModuleCat R) G` を繋ぐ比較関手
4. 局所 Tate 双対性そのもの(`absent-recheck` で 0 件のまま)

## 逸脱の記録

* **逸脱 1**: 原典の `M` は `Zp`-加群として `Z/pnZ` で、`Γ_K` が連続に作用する。
  本ファイルの `charRep ψ` は `k = ZMod (p^n)` 上の**階数 1 の表現**として同じものを作るが、
  **連続性を課していない**(離散群コホモロジーを使う)。
  ★これは (a)(b) を測るための層であり、原典の `H²` そのものではない。
  ★**この差は `discreteCardH2` の名前と docstring に明記した。**
* **逸脱 2**: 原典は `H2(K, M)` の**同型類**を問題にするが、本ファイルの
  `cardCohomologyChar` は**位数**である。`LocalTateDualityRoute` の逸脱 1 と同じ理由
  (位数に落とすと `H²` を定義せずに済む)であり、そちらと整合している。

## 退化の自己検査

* ★`natCard_groupCohomology_zero` + `invariants_charRep_one` ——
  自明指標のとき `H⁰` はちょうど `k` である。**有限性インスタンスは空虚ではない**
  (`Finite` を付けた対象が `0` でないことを実際に見ている)。
* ★`natCard_groupCohomology_le` —— 上界が `|A|^(|G|^n)` であって `0` や `1` ではない。
* ★`groupCohomologyMulEquivIso` は `Iso` であって単なる射ではない
  (`hom_inv_id` / `inv_hom_id` を両方証明した)。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open CategoryTheory groupCohomology

universe u v

/-! ## 1. 抽象核 —— 語彙ゼロ

★この節には表現論もコホモロジーも出てこない。純粋な論理と圏論だけである。 -/

section AbstractCore

/-- **★抽象核 0** —— 「`f` に沿った帰納法原理」から `f` の全射性が出る。

`h` は「`C (f a)` が全ての `a` で成り立てば `C` は恒真」という形をした帰納法原理である。
mathlib はコホモロジーの元について**全射性ではなく帰納法原理**の形で補題を用意している
(`groupCohomology.H2_induction_on`)ので、その形から全射性を取り出す一行が要る。

★分岐・付値・Galois・コホモロジー・表現論の語彙が 1 語も出ない。
`α β` は `Sort` でよく、構造は一切要らない。 -/
theorem surjective_of_inductionOn {α : Sort*} {β : Sort*} (f : α → β)
    (h : ∀ C : β → Prop, (∀ a : α, C (f a)) → ∀ b : β, C b) : Function.Surjective f :=
  h (fun b => ∃ a, f a = b) (fun a => ⟨a, rfl⟩)

/-- **★抽象核 0(系)** —— 有限型からの「帰納法原理」があれば行き先も有限。 -/
theorem finite_of_inductionOn {α β : Type*} [Finite α] (f : α → β)
    (h : ∀ C : β → Prop, (∀ a : α, C (f a)) → ∀ b : β, C b) : Finite β :=
  Finite.of_surjective f (surjective_of_inductionOn f h)

/-- **★抽象核 1** —— `ModuleCat` の同型は台集合の濃度を保つ。

★具体圏の忘却関手が同型を同型に送ることだけを使う。表現論もコホモロジーも出てこない。 -/
theorem natCard_eq_of_moduleCatIso {k : Type u} [CommRing k] {X Y : ModuleCat.{v} k}
    (i : X ≅ Y) : Nat.card X = Nat.card Y :=
  Nat.card_congr (((forget (ModuleCat.{v} k)).mapIso i).toEquiv)

end AbstractCore

/-! ## 2. 抽象核 —— 表現の制限の合成(コホモロジーの語彙ゼロ)

★★ここが (b) の核である。**群同型は要らない**。
`f : G →* H` と `g : H →* G` について `f ∘ g = id` という**片側の等式だけ**で
`res g (res f A) ≅ A` が出る。★原典(および欲しかった形)より一般で、しかも易しい。 -/

section ResCore

variable {k G H : Type} [CommRing k] [Group G] [Group H]

/-- **★★抽象核 2** —— `f.comp g = id` なら `res g (res f A) ≅ A`。

★**群同型は仮定していない**。片側の恒等式だけでよい。
台加群は同じで、作用が `A.ρ (f (g x))` から `A.ρ x` に変わるだけなので、
恒等線形写像がそのまま絡み合い写像になる。 -/
noncomputable def resResIso {f : G →* H} {g : H →* G} (h : f.comp g = MonoidHom.id H)
    (A : Rep k H) : Rep.res g (Rep.res f A) ≅ A :=
  Rep.mkIso (Representation.Equiv.mk (LinearEquiv.refl k A) (fun x => by
    have hx : f (g x) = x := by rw [← MonoidHom.comp_apply, h]; rfl
    simp [hx]))

/-- **★抽象核 3** —— `groupCohomology.map` の引数についての合同則。

mathlib の `groupCohomology.congr` は `h ▸ φ` という輸送を含むので、
そのままでは `φ` の側を比較できない。**`f₁ f₂` を変数のまま量化して `subst` を使える形**に
書き直すと、`φ` の比較は台写像の一致だけで済む。

★`lean-idioms.md` の「一般化してから `subst`」の形。 -/
theorem map_congr {A : Rep k H} {B : Rep k G} {f₁ f₂ : G →* H} (h : f₁ = f₂)
    {φ₁ : Rep.res f₁ A ⟶ B} {φ₂ : Rep.res f₂ A ⟶ B} (hφ : ∀ x : A, φ₁.hom x = φ₂.hom x)
    (n : ℕ) : groupCohomology.map f₁ φ₁ n = groupCohomology.map f₂ φ₂ n := by
  subst h
  have : φ₁ = φ₂ := by ext x; exact hφ x
  rw [this]

end ResCore

/-! ## 3. ★★(a) —— 有限群・有限係数の群コホモロジーは有限

★原典が使うのは `n = 2` だけだが、証明は次数に一切触れないので全次数で通る。

道筋は 3 段:
`(Fin n → G) → A` が有限 ⟹ `cocycles A n`(核なので単射で入る)が有限
⟹ `groupCohomology A n`(`π` が全射)が有限。 -/

section Finiteness

variable {k G : Type} [CommRing k] [Group G]

/-- 非斉次コチェイン `C^n(G, A) = ((Fin n → G) → A)` は有限。 -/
instance finiteCochains (A : Rep k G) (n : ℕ) [Finite G] [Finite A] :
    Finite ((inhomogeneousCochains A).X n) := by
  show Finite ((Fin n → G) → A)
  infer_instance

/-- コサイクル `Z^n(G, A)` は有限(コチェインへ単射で入るから)。 -/
instance finiteCocycles (A : Rep k G) (n : ℕ) [Finite G] [Finite A] :
    Finite (groupCohomology.cocycles A n) :=
  Finite.of_injective (ConcreteCategory.hom (groupCohomology.iCocycles A n))
    ((ModuleCat.mono_iff_injective _).1 inferInstance)

/-- **★★(a) —— 有限群の有限係数コホモロジーは有限**。

★これが `LocalTateDualityRoute.lean` の測定表で
「`Finite (groupCohomology.H2 A)` は `failed to synthesize`」と記録された欠落である。
★`H2` は `groupCohomology A 2` の `abbrev` なので、これで `Finite ↑(H2 A)` も
`inferInstance` で出る。 -/
instance finiteGroupCohomology (A : Rep k G) (n : ℕ) [Finite G] [Finite A] :
    Finite (groupCohomology A n) :=
  Finite.of_surjective (ConcreteCategory.hom (groupCohomology.π A n))
    ((ModuleCat.epi_iff_surjective _).1
      (HomologicalComplex.instEpiHomologyπ (inhomogeneousCochains A) n))

/-- **★(a) の名指しの形** —— 原典が要求する `H²` の有限性。 -/
theorem finite_H2 (A : Rep k G) [Finite G] [Finite A] :
    Finite (groupCohomology.H2 A) :=
  inferInstance

/-- **★`H2π` は全射**。

mathlib は `groupCohomology.H2_induction_on`(帰納法原理)しか持っていない。
抽象核 0 を当てて全射性に直したもの。★低次の明示表示を使う道でも (a) が出る。 -/
theorem H2π_surjective (A : Rep k G) :
    Function.Surjective (ConcreteCategory.hom (groupCohomology.H2π A)) :=
  surjective_of_inductionOn _ (fun _C hC x => groupCohomology.H2_induction_on x hC)

/-- **★濃度の上界** —— `|H^n(G, A)| ≤ |A|^(|G|^n)`。

★退化の自己検査でもある。上界が `0` や `1` ではないことが見える。 -/
theorem natCard_groupCohomology_le (A : Rep k G) (n : ℕ) [Finite G] [Finite A] :
    Nat.card (groupCohomology A n) ≤ Nat.card A ^ Nat.card G ^ n := by
  refine le_trans (Nat.card_le_card_of_surjective _
    ((ModuleCat.epi_iff_surjective _).1
      (HomologicalComplex.instEpiHomologyπ (inhomogeneousCochains A) n))) ?_
  refine le_trans (Nat.card_le_card_of_injective
    (β := ((inhomogeneousCochains A).X n)) _
    ((ModuleCat.mono_iff_injective (groupCohomology.iCocycles A n)).1 inferInstance)) ?_
  have hX : Nat.card ((inhomogeneousCochains A).X n) = Nat.card ((Fin n → G) → A) := rfl
  rw [hX, Nat.card_fun, Nat.card_fun, Nat.card_eq_fintype_card (α := Fin n), Fintype.card_fin]

/-- **★`H⁰` は不変元** —— 退化の自己検査に使う。 -/
theorem natCard_groupCohomology_zero (A : Rep k G) :
    Nat.card (groupCohomology A 0) = Nat.card A.ρ.invariants :=
  natCard_eq_of_moduleCatIso (groupCohomology.H0Iso A)

end Finiteness

/-! ## 4. ★★(b) —— 群同型に沿ったコホモロジーの同型

★★原典の「This is clearly a group-theoretic condition on M」の中身である。

原文 (pGC p.3):
> This is clearly a group-theoretic condition on M.

★原典が要求するのは `n = 2` だけだが、証明は次数に触れないので全次数で通る。 -/

section MulEquivInvariance

variable {k G H : Type} [CommRing k] [Group G] [Group H]

/-- `e : G ≃* H` に沿った制限の合成は恒等(抽象核 2 の特殊化)。 -/
noncomputable def resMulEquivIso (e : G ≃* H) (A : Rep k H) :
    Rep.res e.symm.toMonoidHom (Rep.res e.toMonoidHom A) ≅ A :=
  resResIso (by ext x; simp) A

/-- **★★(b) —— 群同型に沿ったコホモロジーの同型(全次数)**。

★これが `LocalTateDualityRoute.lean` の測定表で
「`e : G ≃* H` について `H2 (Rep.res e.toMonoidHom A) ≅ H2 A` は `exact?` が閉じられない。
`groupCohomology.map` と `congr` から組めるはずだが束ねられていない」と
記録された欠落である。

★2 つの合成が恒等になることを `map_comp` + 抽象核 3(`map_congr`)+ `map_id` で示す。 -/
noncomputable def groupCohomologyMulEquivIso (e : G ≃* H) (A : Rep k H) (n : ℕ) :
    groupCohomology A n ≅ groupCohomology (Rep.res e.toMonoidHom A) n where
  hom := groupCohomology.map e.toMonoidHom (𝟙 (Rep.res e.toMonoidHom A)) n
  inv := groupCohomology.map e.symm.toMonoidHom (resMulEquivIso e A).hom n
  hom_inv_id := by
    rw [← groupCohomology.map_comp]
    exact Eq.trans (map_congr (A := A) (B := A)
      (f₁ := e.toMonoidHom.comp e.symm.toMonoidHom) (f₂ := MonoidHom.id H)
      (φ₁ := (Rep.resFunctor e.symm.toMonoidHom).map (𝟙 (Rep.res e.toMonoidHom A)) ≫
        (resMulEquivIso e A).hom)
      (φ₂ := 𝟙 A) (by ext x; simp) (fun _ => rfl) n)
      (groupCohomology.map_id n)
  inv_hom_id := by
    rw [← groupCohomology.map_comp]
    exact Eq.trans (map_congr (A := Rep.res e.toMonoidHom A) (B := Rep.res e.toMonoidHom A)
      (f₁ := e.symm.toMonoidHom.comp e.toMonoidHom) (f₂ := MonoidHom.id G)
      (φ₁ := (Rep.resFunctor e.toMonoidHom).map (resMulEquivIso e A).hom ≫
        𝟙 (Rep.res e.toMonoidHom A))
      (φ₂ := 𝟙 (Rep.res e.toMonoidHom A)) (by ext x; simp) (fun _ => rfl) n)
      (groupCohomology.map_id n)

/-- **★(b) の名指しの形** —— 原典が要求する `H²` の同型。 -/
noncomputable def H2MulEquivIso (e : G ≃* H) (A : Rep k H) :
    groupCohomology.H2 A ≅ groupCohomology.H2 (Rep.res e.toMonoidHom A) :=
  groupCohomologyMulEquivIso e A 2

/-- **★(b) を位数で述べたもの** —— 消費側(`LocalTateDualityData`)が要るのはこの形。 -/
theorem natCard_groupCohomology_eq_of_mulEquiv (e : G ≃* H) (A : Rep k H) (n : ℕ) :
    Nat.card (groupCohomology (Rep.res e.toMonoidHom A) n) = Nat.card (groupCohomology A n) :=
  (natCard_eq_of_moduleCatIso (groupCohomologyMulEquivIso e A n)).symm

end MulEquivInvariance

/-! ## 5. 指標に付随する階数 1 の表現

原典の `M`(`Zp`-加群として `Z/pnZ`、`Γ_K` は指標で作用)にあたる層。
★一般の可換環 `k` と一般の群 `G` で書く。`p` も `Γ_K` も出てこない。 -/

section CharRep

variable {k G H : Type} [CommRing k] [Group G] [Group H]

/-- **指標 `ψ : G →* kˣ` に付随する階数 1 の表現**。

台加群は `k` 自身で、`g` は `ψ g` 倍として作用する。

★原典の `M`(「`M` is isomorphic as a `Zp`-module to `Z/pnZ`」)の一般形である。
`k = ZMod (p^n)` と取れば原典の設定になる。 -/
noncomputable def charRep (ψ : G →* kˣ) : Rep k G :=
  Rep.of ((Algebra.lmul k k).toRingHom.toMonoidHom.comp ((Units.coeHom k).comp ψ))

@[simp] theorem charRep_ρ_apply (ψ : G →* kˣ) (g : G) (x : k) :
    (charRep ψ).ρ g x = (ψ g : k) * x := rfl

/-- **★制限は指標の合成**。★`rfl` で通る(定義通り)。 -/
theorem res_charRep (f : G →* H) (ψ : H →* kˣ) :
    Rep.res f (charRep ψ) = charRep (ψ.comp f) := rfl

instance finiteCharRep (ψ : G →* kˣ) [Finite k] : Finite (charRep ψ) := by
  show Finite k
  infer_instance

/-- 指標 `ψ` に付随する `H^n` の位数。 -/
noncomputable def cardCohomologyChar (n : ℕ) (ψ : G →* kˣ) : ℕ :=
  Nat.card (groupCohomology (charRep ψ) n)

/-- **★★★原典の「群論性」は定理である**。

原文 (pGC p.3):
> This is clearly a group-theoretic condition on M.

★群同型 `e : G ≃* H` に沿って指標を引き戻すと `H^n` の位数は変わらない。
★(b)(`groupCohomologyMulEquivIso`)と `res_charRep`(`rfl`)だけから出る。
★有限性は要らない(`Nat.card` は無限でも定義されている)。 -/
theorem cardCohomologyChar_comp_mulEquiv (e : G ≃* H) (ψ : H →* kˣ) (n : ℕ) :
    cardCohomologyChar n (ψ.comp e.toMonoidHom) = cardCohomologyChar n ψ :=
  natCard_groupCohomology_eq_of_mulEquiv e (charRep ψ) n

/-- **★退化の自己検査** —— 自明指標の不変元は全体。 -/
theorem invariants_charRep_one : (charRep (1 : G →* kˣ)).ρ.invariants = ⊤ :=
  Submodule.eq_top_iff'.2 (fun x => (Representation.mem_invariants _ x).2 (fun _ => by simp))

/-- **★★壁は空虚ではない** —— 自明指標のとき `H⁰` はちょうど `k` である。

★(a) で作った有限性インスタンスが `0` でない対象に付いていることの確認。 -/
theorem cardCohomologyChar_zero_one :
    cardCohomologyChar 0 (1 : G →* kˣ) = Nat.card k := by
  rw [cardCohomologyChar, natCard_groupCohomology_zero, invariants_charRep_one]
  exact Nat.card_congr (Equiv.subtypeUnivEquiv (fun _ => trivial))

end CharRep

/-! ## 6. ★★(c) の到達点 —— inflation は全次数で在る。無いのは完全性と colimit

★★**2026-09-04 の記録「inflation-restriction は次数 1 だけ」は不正確だった。**
`groupCohomology.infNatTrans` は**全次数の自然変換**であり、`n = 2` でも
`H²(G/S, A^S) ⟶ H²(G, A)` という射が実際に作れる(下の `inflation`)。

★**無いのは次の 2 つである**:
1. `n ≥ 2` での完全性(`H1InfRes_exact` は `n = 1` だけ。`H2InfRes` は `Unknown constant`)
2. 開正規部分群の有向系についての `colim_S H^n(G/S, A^S) ≅ H^n_cont(G, A)`
   (`groupCohomology.colimitIso` は `Unknown constant`)

★したがって「有限次で得た情報を `Γ_K` に上げる」道は、
**射は引けるが、それが同型であることが言えない**という段階にある。 -/

section Inflation

variable {k G : Type} [CommRing k] [Group G]

/-- **★inflation 写像(全次数)**。

`S ⊴ G` について `H^n(G/S, A^S) ⟶ H^n(G, A)`。
mathlib の `groupCohomology.infNatTrans` を対象で評価しただけのものだが、
「`H²` の inflation が無い」という 2026-09-04 の記録を訂正する測定結果なので名前を付けておく。

★★**これが同型であること(あるいは完全列に入ること)は `n = 1` でしか証明されていない。**
`n = 2` の完全性が (c) の 1 つ目の欠落である。 -/
noncomputable def inflation (A : Rep k G) (S : Subgroup G) [S.Normal] (n : ℕ) :
    groupCohomology ((Rep.quotientToInvariantsFunctor k S).obj A) n ⟶ groupCohomology A n :=
  (groupCohomology.infNatTrans k S n).app A

/-- **★`n = 2` でも inflation は引ける**(2026-09-07 の測定を宣言として残したもの)。 -/
noncomputable def inflationH2 (A : Rep k G) (S : Subgroup G) [S.Normal] :
    groupCohomology.H2 ((Rep.quotientToInvariantsFunctor k S).obj A) ⟶ groupCohomology.H2 A :=
  inflation A S 2

end Inflation

/-! ## 7. ★消費側(`LocalTateDualityData`)との接続

★★**ここは `LocalTateDualityData` を作る節ではない。**
`isGroupTheoretic` フィールドと**同じ形の命題**が、離散群コホモロジーについては
**定理である**ことを示す節である。 -/

section Consumer

variable {p : ℕ} [Fact p.Prime]

/-- **離散群コホモロジーで測った `|H²(Γ_K, M_ψ)|`**。

★★**これは原典の `H²` では「ない」。** 原典の `H²` は `Γ_K` の位相を見る
連続コホモロジーであり、本定義は位相を無視した抽象群コホモロジーである。
mathlib に副有限群の連続 `H²` が無い(上の (c) の測定)ので、
`LocalTateDualityData.cardH2` の**形だけ**を満たす候補として置いている。

★この候補は `isGroupTheoretic` を満たす(下の定理)が、
`cardH2_eq_natCard`(局所 Tate 双対性)は満たさない。 -/
noncomputable def discreteCardH2 (_K : PAdicLocalField p) (n : ℕ)
    (ψ : _K.absGal →* (ZMod (p ^ n))ˣ) : ℕ :=
  cardCohomologyChar 2 ψ

def discreteCardH2_isGroupTheoretic.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }

/-- **★★★原典の「群論性」の主張を、消費側と同じ形で証明したもの**。

原文 (pGC p.3):
> This is clearly a group-theoretic condition on M.

> Thus, we conclude that the isomorphism class of the ΓK-module Zp(1) can be recovered
> group-theoretically from ΓK.

★これは `LocalTateDualityData.isGroupTheoretic` フィールドと**同じ形**である
(そちらは核が開という仮定を持つが、本定理は**それすら要らない**ので強い)。

★★**ただし本定理は `LocalTateDualityData` を作らない。**
`discreteCardH2` は連続コホモロジーではないので `cardH2_eq_natCard` を満たさないからである。
本定理が言うのは「3 つのフィールドのうち `isGroupTheoretic` は
**H² が関手的に定義されている限り自動である**」ということであり、
★**残る障害はすべて (c)(連続コホモロジーと双対性)に集約される**、という測定結果である。 -/
theorem discreteCardH2_isGroupTheoretic (K K' : PAdicLocalField p)
    (α : ContinuousMulEquiv K.absGal K'.absGal) (n : ℕ)
    (ψ : K'.absGal →* (ZMod (p ^ n))ˣ) :
    discreteCardH2 K n (ψ.comp (α : K.absGal ≃* K'.absGal).toMonoidHom)
      = discreteCardH2 K' n ψ :=
  cardCohomologyChar_comp_mulEquiv _ ψ 2

end Consumer

end ABC3.Found.PGC
