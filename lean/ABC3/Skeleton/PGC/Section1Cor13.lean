import ABC3.Skeleton.PGC.Setup
import ABC3.Interface.PGC.LocalFieldData

/-!
# [pGC] Corollary 1.3 — 惰性群

## 設計の経緯(必読)

惰性群 `I_K` を **`Interface` の自由なデータ**として置く素朴な設計は**壊れている**。
`ABC3/Check/PGC/InertiaDegeneracy.lean` で機械的に示したとおり、
`I_K := ⊥` でも `I_K := ⊤` でも Corollary 1.3 の形は成り立ってしまう
(どちらも `sorry` 無しで証明済み)。「偽になる」のではなく「**自明になる**」型の失敗。

したがって `I_K` は **posit するものではなく derive するもの**である。
原文 p.3 の不分岐判定から構成する:

```
I_K = ⋂ { H ⊆ Γ_K : H は開、かつ q_{L_H} = q^[Γ_K:H] }
```

**この構成自体は原文にない**(暗黙の段 #7)。原文が与えるのは判定条件までで、
そこから `I_K` を得る操作は書かれていない。

## 退化は消えていない——移動した

`I_K` の自由度は消えたが、`RD : ResidueCardinality` の自由度は残っている。
実際、`RD.card := fun _ => p`(`isPrimePow` は f = 1 で満たす)と置くと
不分岐条件は `p = p^[Γ_K:H]` となり、成り立つのは指数 1 のときだけ、
すなわち `I_K = ⊤` に退化する。

これは設計の失敗ではなく**設計どおり**である——退化の可能性が
**1 箇所(`Interface`)に局在した**。Track B が本物の `ResidueCardinality` を
構成した時点で、ここに依存する全ての statement が一斉に非空虚性の検査を受ける。
自由なデータのままなら、退化は各 statement に散らばったまま検査されなかった。
-/

namespace ABC3.Skeleton.PGC

open ABC3.Meta ABC3.Interface.PGC

variable {p : ℕ} [Fact p.Prime]

/-- 開部分群 H に対応する拡大 L_H/K が**不分岐**であることの判定条件。

原文 (pGC p.3):
> Note, moreover, that L is unramified over K if and only q_L = q^[Γ_K:H].

原文の "if and only" は "if and only if" の脱字と読めるが、逐語欄では原文どおりに保つ。 -/
def IsUnramifiedAt (RD : ResidueCardinality p) (SC : SubgroupCorrespondence p)
    (K : PAdicLocalField p) (H : Subgroup K.absGal) (hH : IsOpen (H : Set K.absGal)) : Prop :=
  RD.card (SC.field K H hH) = RD.card K ^ H.index

def IsUnramifiedAt.src : Source :=
  { paper := "pGC", pdfPage := 3, item := "Section 1 (unramifiedness criterion)",
    sectionId := "setup-1-unram" }

/-- 惰性群 `I_K`——不分岐な開部分群すべての共通部分として**構成する**。

★**原文に無い段**(暗黙の段 #7)。原文は判定条件までしか与えず、
そこから `I_K` を得る操作を書いていない。ここは我々の構成であり、
原典の転写ではない。 -/
noncomputable def inertia (RD : ResidueCardinality p) (SC : SubgroupCorrespondence p)
    (K : PAdicLocalField p) : Subgroup K.absGal :=
  sInf {H | ∃ hH : IsOpen (H : Set K.absGal), IsUnramifiedAt RD SC K H hH}

def inertia.src : Source :=
  { paper := "pGC", pdfPage := 3, item := "Corollary 1.3", sectionId := "cor-1-3" }

/-- 惰性群を「K に付随する対象」として束ねたもの。移送は部分群の押し出し。

`inertiaObjectOf`(`Check/PGC/InertiaDegeneracy.lean`)と同じ形だが、
`obj` が**自由なデータでなく構成**である点が本質的な違い。 -/
noncomputable def inertiaObject (RD : ResidueCardinality p) (SC : SubgroupCorrespondence p) :
    AssociatedObject p where
  Obj := fun K => Subgroup K.absGal
  obj := fun K => inertia RD SC K
  transport := fun α S => S.map α.toMulEquiv.toMonoidHom

def inertiaObject.src : Source :=
  { paper := "pGC", pdfPage := 3, item := "Corollary 1.3", sectionId := "cor-1-3" }

/-- **[pGC] Corollary 1.3**

原文 (pGC p.3):
> The inertia subgroup I_K ⊆ Γ_K can be determined group-theoretically from Γ_K.

## 条件付き形式化

`RD`(剰余体の元の個数)と `SC`(開部分群と中間体の対応)を仮説に取る。
どちらも原文が §1 で所与として使うが、我々がまだ構成できていないもの。

## 未解決

原文の論拠は Proposition 1.2 を (L, H) に適用することだが、
そこから `I_K` そのものを得る段は書かれていない(上記の設計の経緯を参照)。 -/
theorem inertia_recoverable (RD : ResidueCardinality p) (SC : SubgroupCorrespondence p) :
    (inertiaObject RD SC).RecoverableFromAbsGal := by
  sorry

def inertia_recoverable.src : Source :=
  { paper := "pGC", pdfPage := 3, item := "Corollary 1.3", sectionId := "cor-1-3" }

/-- 原文の証明文から抽出した、証明が要求するもの(G6)。★下界。 -/
def inertia_recoverable.needs : List ProofObligation :=
  [ .derivation "Proposition 1.2 を (L, H) に適用して q_L を得る段" 3,
    .implicitStep
      "不分岐判定から惰性群 I_K そのものを構成する段(原文に無い。ここでは開かつ不分岐な部分群の共通部分として構成した)" 3 ]

end ABC3.Skeleton.PGC
