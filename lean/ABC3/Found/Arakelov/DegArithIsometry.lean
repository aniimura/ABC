/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.DegArithPre
import ABC3.Meta.Claim

/-!
# `deg_F` は**等長同型で不変**である（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

## ★★★★★★★★これは何か

`§9-780` で `degArithPre` は**切断の取り方に依らず**（段 D）、**加法的**（段 E）になった。
本ファイルは残る一点——**等長同型で不変**であること——を入れる。
★これで `deg_F` が `APicM`（等長同型類の群）の上の関数になる用意が整う。

## ★★★★★機構（2 段）

    norm_of_isIsometry     : 等長同型は切断のノルムを保つ
    degFinPre_of_iso       : 有限部分は前層の同型で不変（`Γ(L) ≅ Γ(M)` を `⊤` で取るだけ）

★`IsIsometry L M φ` の定義は「`L.h V (pullTriv φ V e) p = M.h V e p`」であり、
在庫の `trivValue_pullTriv`（`trivValue L V (pullTriv φ V e) s = trivValue M V e (φ s)`）と
合わせると**切断のノルムがそのまま等しい**ことが出る。

★★有限部分は `φ` を `op ⊤` で評価した `Γ(X,⊤)`-線形同型が
`Γ(X,⊤)·s` を `Γ(X,⊤)·(φ s)` へ移すことに尽きる。

## ★残っている段（明示）

★★`deg_F` を `APicM` の上の関数として**定義する**には、
さらに「非零切断のノルムはどの複素点でも非零」（`degArithPre` の欄 `hnorm`）が要る。
★★★`Γ(L)` は可逆なので `s ≠ 0` ならイデアルの非零元に対応し、
埋め込みの像も非零のはずであるが、その段はまだ入っていない。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite NumberField
open ABC3.Found.GenEll

/-- ★★★★★★★★**等長同型は切断のノルムを保つ**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★`IsIsometry` の定義（基準ノルム `h` が `pullTriv` で対応する）と、
在庫の `trivValue_pullTriv`（自明化での値が `φ` と可換）を合わせるだけである。 -/
theorem norm_of_isIsometry {X : Scheme.{0}} {L M : AMetric X} {φ : L.sheaf ≅ M.sheaf}
    (hφ : IsIsometry L M φ) (s : (L.sheaf.obj (op ⊤) : Type))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) :
    M.norm ((φ.hom.app (op ⊤)) s) p = L.norm s p := by
  obtain ⟨c⟩ := nonempty_normChart M.triv p
  rw [AMetric.norm_eq M _ p c.V c.e c.hp,
      AMetric.norm_eq L s p c.V (pullTriv φ c.V c.e) c.hp]
  show trivSecNorm M.sheaf c.V c.e p c.hp _ * M.metric.h c.V c.e p
    = trivSecNorm L.sheaf c.V (pullTriv φ c.V c.e) p c.hp s
      * L.metric.h c.V (pullTriv φ c.V c.e) p
  rw [hφ c.V c.e p c.hp]
  congr 1
  show ‖evalOn p c.V c.hp (trivValue M.sheaf c.V c.e _)‖
    = ‖evalOn p c.V c.hp (trivValue L.sheaf c.V (pullTriv φ c.V c.e) s)‖
  rw [trivValue_pullTriv]

/-- ★★★★★**有限部分は前層の同型で不変**である。

★`φ` を `op ⊤` で評価した `Γ(X,⊤)`-線形同型が `Γ(X,⊤)·s` を `Γ(X,⊤)·(φ s)` へ移す。 -/
theorem degFinPre_of_iso {X : Scheme.{0}} (L M : AInv X)
    (φ : L.carrier.sheaf ≅ M.carrier.sheaf) (s : (L.carrier.sheaf.obj (op ⊤) : Type)) :
    degFinPre M ((φ.hom.app (op ⊤)) s) = degFinPre L s := by
  set e : (L.carrier.sheaf.obj (op ⊤) : Type) ≃ₗ[(Γ(X, (⊤ : X.Opens)) : Type)]
      (M.carrier.sheaf.obj (op ⊤) : Type) :=
    ((PresheafOfModules.evaluation X.ringCatSheaf.obj (op ⊤)).mapIso φ).toLinearEquiv with he
  have hmap : Submodule.map (e : _ →ₗ[(Γ(X, (⊤ : X.Opens)) : Type)] _)
      ((Γ(X, (⊤ : X.Opens)) : Type) ∙ s)
      = (Γ(X, (⊤ : X.Opens)) : Type) ∙ ((φ.hom.app (op ⊤)) s) := by
    rw [Submodule.map_span]
    simp only [Set.image_singleton, LinearEquiv.coe_coe]
    rfl
  show Real.log (Nat.card (_ ⧸ ((Γ(X, (⊤ : X.Opens)) : Type) ∙ ((φ.hom.app (op ⊤)) s)))) = _
  congr 2
  exact (Nat.card_congr (Submodule.Quotient.equiv _ _ e hmap).toEquiv).symm

/-- ★★★★★★★★★★**`deg_F` は等長同型で不変である**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★★これで `deg_F` が `APicM`（等長同型類の群）の上の関数になる用意が整った。 -/
theorem degArithPre_of_isIsometry (F : Type) [Field F] [NumberField F]
    (L M : AInv (Spec (CommRingCat.of (𝓞 F))))
    {φ : L.carrier.sheaf ≅ M.carrier.sheaf} (hφ : IsIsometry L.carrier M.carrier φ)
    (s : (L.carrier.sheaf.obj (op ⊤) : Type)) :
    degArithPre F M ((φ.hom.app (op ⊤)) s) = degArithPre F L s := by
  have hA : archDeg F M.carrier ((φ.hom.app (op ⊤)) s) = archDeg F L.carrier s := by
    show -(∑ σ : F →+* ℂ,
        Real.log (M.carrier.norm ((φ.hom.app (op ⊤)) s) (embSpecPoint F σ))) / _
      = -(∑ σ : F →+* ℂ, Real.log (L.carrier.norm s (embSpecPoint F σ))) / _
    simp only [norm_of_isIsometry hφ]
  show degFinPre M _ / _ + archDeg F M.carrier _ = _
  rw [degFinPre_of_iso L M φ s, hA]
  rfl

/-! ### ★出典の紐付け(`.src`) -/

def norm_of_isIsometry.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(等長同型は切断のノルムを保つ)",
    sectionId := "genell-def-1-1-ii" }

def degArithPre_of_isIsometry.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(deg_F は等長同型で不変)",
    sectionId := "genell-def-1-1-ii" }

def degArithPre_of_isIsometry.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "trivValue_pullTriv(自明化での値は φ と可換、在庫)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.trivValue_pullTriv") 3,
    .citation "[ABC3]" "AMetric.norm_eq(ノルムはチャートの取り方に依らない、在庫)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.AMetric.norm_eq") 3,
    .implicitStep
      ("★★deg_F を APicM の上の関数として**定義する**には、" ++
       "さらに「非零切断のノルムはどの複素点でも非零」が要る。" ++
       "★★★Γ(L) は可逆なので s ≠ 0 ならイデアルの非零元に対応し、" ++
       "埋め込みの像も非零のはずであるが、その段はまだ入っていない") 4 ]

end ABC3.Found.Arakelov
