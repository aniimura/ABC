import ABC3.Found.GenEll.HeightNonneg
import ABC3.Found.GenEll.ProductFormula

/-!
# [GenEll] Proposition 1.4, (iii) —— 高さは**類にしか依らない**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

## ★★★(iii) の機構は「次数が `APic` に降りる」ことである

原文 (iii) は

> (iii) The BD-class of htL depends only on the isomorphism class of the

すなわち「高さの BD-class は `L_ℚ` の同型類だけに依る」である。

★★構成の側で見ると、これは **`deg` が主算術因子の上で消える**(積公式)ことの
帰結にすぎない。`Spec 𝓞_F` 上では

    ADiv(F)/APrc(F) ≅ APic(Spec 𝓞_F)

であり、`degNormalizedAPic` が既に構成してある(`ProductFormula.lean`)。
★★★**したがって `htArith` は `APic(Spec 𝓞_F)` を経由する。**

## ★★★原文より強い形が出る

原文の結論は **BD-class の一致**(定数差を許す)である。
★★構成の側では **`ht` そのものが一致する**——定数差すら出ない。
★理由: 高さは `x_F^*` の**類**だけで決まり、類が同じなら次数は同じだからである。

## ★ここにも `PullbackMul` は要らない

(i)(加法性)は引き戻しの積の保存を要したが、
**(iii) は因子を 1 つずつ見るだけなので要らない**。★無条件である。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField

variable (F : Type) [Field F] [NumberField F]

/-! ## ★★高さは `APic` を経由する -/

/-- ★★**引き戻しの算術因子類** —— `x_F^* D̄` の `APic(Spec 𝓞_F)` における類。

★原文の `APic(Spec(O_F))` の元そのものである。 -/
noncomputable def pullbackAPic {X : Scheme.{0}} (D : ArithCartier X)
    (xF : specRingOfIntegers F ⟶ X) : APicOF F :=
  QuotientAddGroup.mk' (APrc F) (pullbackADiv F D xF)

/-- ★★★**高さは類の次数である** —— `ht = deg_APic ∘ (類を取る)`。

★★これが「高さが同型類にしか依らない」ことの内容である。 -/
theorem htArith_eq_degNormalizedAPic {X : Scheme.{0}} (D : ArithCartier X)
    (xF : specRingOfIntegers F ⟶ X) :
    htArith F D xF = degNormalizedAPic (pullbackAPic F D xF) :=
  (degNormalizedAPic_mk _).symm

/-! ## ★★★`Proposition 1.4, (iii)` -/

/-- ★★★**類が同じなら高さは等しい**。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

★★★**原文は BD-class の一致(定数差を許す)を言うが、
構成の側では `ht` そのものが一致する。** 定数差すら出ない。

★機構は `deg` が主算術因子の上で消えること(積公式、`ProductFormula.lean`)。 -/
theorem htArith_congr_of_pullbackAPic_eq {X : Scheme.{0}} (D E : ArithCartier X)
    (xF : specRingOfIntegers F ⟶ X)
    (h : pullbackAPic F D xF = pullbackAPic F E xF) :
    htArith F D xF = htArith F E xF := by
  rw [htArith_eq_degNormalizedAPic, htArith_eq_degNormalizedAPic, h]

/-- ★★**主算術因子だけ違うなら高さは等しい** —— 上の命題の具体形。

★`x_F^*D − x_F^*E` が主算術因子(= ある `f ∈ Fˣ` の因子)なら高さは一致する。 -/
theorem htArith_congr_of_sub_mem_APrc {X : Scheme.{0}} (D E : ArithCartier X)
    (xF : specRingOfIntegers F ⟶ X)
    (h : pullbackADiv F D xF - pullbackADiv F E xF ∈ APrc F) :
    htArith F D xF = htArith F E xF := by
  refine htArith_congr_of_pullbackAPic_eq F D E xF ?_
  have hz : (QuotientAddGroup.mk' (APrc F))
      (pullbackADiv F D xF - pullbackADiv F E xF) = 0 :=
    (QuotientAddGroup.eq_zero_iff _).2 h
  rw [map_sub, sub_eq_zero] at hz
  exact hz

/-- ★**BD-class の水準でも一致**(原文の結論の形)。

★★定数 `C = 0` で足りる——上で等式そのものが出ているからである。 -/
theorem htArith_bdeq_of_pullbackAPic_eq {X : Scheme.{0}} (D E : ArithCartier X)
    (h : ∀ xF : specRingOfIntegers F ⟶ X, pullbackAPic F D xF = pullbackAPic F E xF) :
    BDeq (fun xF => htArith F D xF) (fun xF => htArith F E xF) :=
  ⟨0, fun xF => by
    show |htArith F D xF - htArith F E xF| ≤ 0
    rw [htArith_congr_of_pullbackAPic_eq F D E xF (h xF), sub_self, abs_zero]⟩

/-! ## ★出典の紐付け(`.src`)

★条つきである。`Proposition 1.4` 全体には (i)〜(iv) の 4 条すべてが要る。 -/

def htArith_congr_of_pullbackAPic_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iii)(構成した高さが APic の類にしか依らないこと)",
    sectionId := "genell-prop-1-4" }

end ABC3.Found.GenEll
